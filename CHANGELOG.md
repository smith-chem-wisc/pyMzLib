# Changelog

Notable changes to pyMzLib. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow semantic versioning judged on the **Python** API — a change to the internal JSON
envelope is not a breaking change unless Python callers can see it.

## [Unreleased]

### Added
- PRIDE Archive support: `pymzlib.pride.list_files()`, `download()`, `total_size_bytes()`,
  and the `PrideFile` type.
- Peptidoforms: `pymzlib.peptidoform.fragments()` — digest an annotated UniProt protein, apply its
  modifications, and fragment every peptide, disclosing the rules it applied (modification census,
  the silent isoform cap, the fixed-charge convention).
- FlashLFQ label-free quantification: `pymzlib.flashlfq.quantify()` — quantify a search's peptides
  across mzML runs with match-between-runs, returning typed `FlashLfqResults` / `Peptide` /
  `ProteinGroup` / `Peak`. Match-between-runs transfers are exposed via `result.peaks`.
- A self-contained .NET payload bundled in the wheel, so no .NET installation is required.
- Documentation site, including the reasoning behind each design decision.

[Unreleased]: https://github.com/smith-chem-wisc/pyMzLib/commits/main