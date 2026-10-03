# quant: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibflashlfq). For worked examples, see the
[FlashLFQ guide](../guides/flashlfq.md).

**Python's keywords are mzLib's parameter names.** The wire abbreviates FlashLFQ's settings
(`ppm`, `mbr`, `mbr-q`, `threads`); `pymzlib.flashlfq.quantify()` spells them as mzLib's
`FlashLfqEngine` does (`ppm_tolerance`, `match_between_runs`, `mbr_q_value_threshold`,
`max_threads`), so a call reads like the C# it runs. The runs travel on stdin, rendered from
`spectra=`; for median polish, the design is rendered from `design=`.

**`median_polish()` returns the protein list only.** The wire's envelope also carries the peptide
file, the samples in fit order, the peptide count and the output directory; Python does not
project those.

## `quantify`

--8<-- "docs/reference/_generated/quant.flashlfq.md"

## `median_polish`

--8<-- "docs/reference/_generated/quant.median-polish.md"
