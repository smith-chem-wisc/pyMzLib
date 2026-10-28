"""Cut a small, real subset of Ensembl release 116's Compara inputs for ``genes map-residues``.

    python make_inputs.py <folder holding Ensembl's release-116 files> <output folder>

The source folder holds the files exactly as Ensembl publishes them:

- ``Compara.116.protein_default.aa.fasta.gz`` (emf/ensembl-compara/homologies/; 908 MB);
- ``vertebrates.GeneTree_content.default.e116.txt.gz`` (compara/);
- ``Homo_sapiens.GRCh38.116.gtf.gz``, ``Mus_musculus.GRCm39.116.gtf.gz`` and
  ``Rattus_norvegicus.GRCr8.116.gtf.gz`` (gtf/<species>/).

Nothing is edited: every line written is a line of the source, byte for byte. What is cut:

- each GTF: its ``#!`` header lines, then the ``gene`` rows of the genes below;
- the gene-tree file: every row (canonical or not) of the genes below, in file order;
- the alignment: the two alignments (gene trees) that hold the canonical protein of human TP53 or
  human GAPDH, in file order. From each, the rows of the genes' canonical proteins are kept, plus
  the first zebrafish row, which the reader checks but does not keep (no gene set holds it). A row
  keeps its full width and every line of it, so each column number is the release's own column.
  Other species' rows, and the blank lines between alignments, are left out; that removes rows,
  never columns.

gzip with mtime 0, no file name in the header and level 9, so two runs write the same bytes.
"""
import gzip
import sys
from pathlib import Path

GENES = {
    "homo_sapiens": ("Homo_sapiens.GRCh38.116.gtf", ["ENSG00000141510", "ENSG00000111640"]),  # TP53, GAPDH
    "mus_musculus": ("Mus_musculus.GRCm39.116.gtf", ["ENSMUSG00000059552", "ENSMUSG00000057666"]),  # Trp53, Gapdh
    "rattus_norvegicus": ("Rattus_norvegicus.GRCr8.116.gtf", ["ENSRNOG00000010756"]),  # Tp53
}
TREES = "vertebrates.GeneTree_content.default.e116.txt.gz"
ALIGNMENT = "Compara.116.protein_default.aa.fasta.gz"
#: The alignments to keep: those holding the canonical protein of human TP53 or GAPDH.
ANCHORS = {"ENSP00000269305", "ENSP00000229239"}
NOT_KEPT_PREFIX = "ENSDARP"  # zebrafish: in the alignment, in no gene set


def gz_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as f:
            f.write(data)


def main() -> None:
    source, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    wanted_genes = {g for _, genes in GENES.values() for g in genes}

    # GTFs: header lines, then the genes' gene rows.
    for name, genes in GENES.values():
        keep = []
        with gzip.open(source / (name + ".gz"), "rb") as f:
            for line in f:
                if line.startswith(b"#!"):
                    keep.append(line)
                    continue
                cells = line.split(b"\t")
                if len(cells) >= 9 and cells[2] == b"gene" and any(
                        f'gene_id "{g}"'.encode() in cells[8] for g in genes):
                    keep.append(line)
        (out / name).write_bytes(b"".join(keep))

    # Gene trees: every row of the genes, and the canonical protein of each.
    tree_rows, canonical = [], set()
    with gzip.open(source / TREES, "rb") as f:
        for line in f:
            cells = line.rstrip(b"\r\n").split(b"\t")
            if len(cells) == 4 and cells[2].decode() in wanted_genes:
                tree_rows.append(line)
                if cells[3] == b"Y":
                    canonical.add(cells[1].decode())
    gz_write(out / TREES, b"".join(tree_rows))

    # Alignment: the rows of those canonical proteins (and one zebrafish row) of the anchors' alignments.
    blocks, rows, row, has_anchor, fish_taken = [], [], None, False, False

    def end_row() -> None:
        nonlocal row, fish_taken
        if row is not None:
            protein = row[0][1:].split()[0].decode()
            if protein in canonical:
                rows.append(row)
            elif protein.startswith(NOT_KEPT_PREFIX) and not fish_taken:
                rows.append(row)
                fish_taken = True
        row = None

    with gzip.open(source / ALIGNMENT, "rb") as f:
        for line in f:
            text = line.strip()
            if text == b"//":
                end_row()
                if has_anchor:
                    blocks.append(b"".join(b"".join(r) for r in rows) + line)
                    if len(blocks) == len(ANCHORS):
                        break
                rows, has_anchor, fish_taken = [], False, False
                continue
            if text.startswith(b">"):
                end_row()
                row = [line]
                if text[1:].split()[0].decode() in ANCHORS:
                    has_anchor = True
                continue
            if row is not None:
                row.append(line)
    gz_write(out / ALIGNMENT, b"".join(blocks))
    print("wrote", sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()))


if __name__ == "__main__":
    main()
