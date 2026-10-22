"""Write mzLib's TestOrthologySnapshot three-species world (mzLib d10017ce) as input files.

Every row is copied from mzLib/Test/FileReadingTests/TestOrthologySnapshot.cs at d10017ce: the
Row() and GeneSet() templates, HumanDump/MouseDump/RatDump and Trees(). Only the file layout
differs: each dump is placed as Ensembl publishes it, homologies/<species>/Compara.116.
protein_default.homologies.tsv.gz (the layout mzLib's ComparaHomologyDump summary names).
"""
import gzip
import sys
from pathlib import Path

out = Path(sys.argv[1])

HUMAN, MOUSE, RAT = "homo_sapiens", "mus_musculus", "rattus_norvegicus"
H1, H2, H3, H4 = "ENSG00000000001", "ENSG00000000002", "ENSG00000000003", "ENSG00000000004"
M1, M2, M3 = "ENSMUSG00000000001", "ENSMUSG00000000002", "ENSMUSG00000000003"
R1, R2 = "ENSRNOG00000000001", "ENSRNOG00000000002"

COLUMNS = [
    "gene_stable_id", "protein_stable_id", "species", "identity", "homology_type",
    "homology_gene_stable_id", "homology_protein_stable_id", "homology_species", "homology_identity",
    "dn", "ds", "goc_score", "wga_coverage", "is_high_confidence", "homology_id",
]
HEADER = "\t".join(COLUMNS) + "\n"


def row(id_, type_, sp, gene, identity, other_sp, other, other_identity, goc="100", high="1"):
    return (f"{gene}\t{gene}P\t{sp}\t{identity}\t{type_}\t{other}\t{other}P\t{other_sp}\t"
            f"{other_identity}\tNULL\tNULL\t{goc}\tNULL\t{high}\t{id_}\n")


def gtf(name, *genes):
    text = "".join(
        f'1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id "{g}"; gene_version "1"; gene_biotype "protein_coding";\n'
        for g in genes)
    (out / name).write_bytes(text.encode("utf-8"))


def gz(relative, text):
    path = out / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 and no file name in the header, so the bytes are the same on every run.
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as f:
            f.write(text.encode("utf-8"))


out.mkdir(parents=True, exist_ok=True)
gtf("Homo_sapiens.GRCh38.116.gtf", H1, H2, H3, H4)
gtf("Mus_musculus.GRCm39.116.gtf", M1, M2, M3)
gtf("Rattus_norvegicus.GRCr8.116.gtf", R1, R2)

dump = "Compara.116.protein_default.homologies.tsv.gz"
gz(f"homologies/{HUMAN}/{dump}", HEADER + row("10", "ortholog_one2one", HUMAN, H1, 80, RAT, R1, 81))
gz(f"homologies/{MOUSE}/{dump}", HEADER
   + row("11", "ortholog_one2one", MOUSE, M1, 90, HUMAN, H1, 91)
   + row("12", "ortholog_one2one", MOUSE, M1, 95, RAT, R1, 96)
   + row("13", "ortholog_one2many", MOUSE, M2, 70, HUMAN, H2, 71, goc="NULL", high="NULL")
   + row("14", "ortholog_one2many", MOUSE, M3, 60, HUMAN, H2, 61, high="0")
   + row("15", "within_species_paralog", MOUSE, M2, 50, MOUSE, M3, 51))
gz(f"homologies/{RAT}/{dump}", HEADER + row("16", "ortholog_one2one", RAT, R2, 88, MOUSE, M2, 89))

gz("vertebrates.GeneTree_content.default.e116.txt.gz",
   f"ENSGT1\t{H1}P\t{H1}\tY\nENSGT1\t{M1}P\t{M1}\tY\nENSGT1\t{R1}P\t{R1}\tY\nENSGT1\t{R1}Q\t{R1}\tN\n"
   f"ENSGT2\t{H2}P\t{H2}\tY\nENSGT2\t{M2}P\t{M2}\tY\nENSGT2\t{M3}P\t{M3}\tY\nENSGT2\t{R2}P\t{R2}\tY\n"
   f"ENSGT3\t{H3}P\t{H3}\tY\nENSGT9\tENSDARP1\tENSDARG00000000001\tY\n")
print("wrote", sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()))
