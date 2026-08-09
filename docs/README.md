# nf-core/starlight: Documentation

The nf-core/starlight documentation is split into the following pages:

- [Architecture](architecture.md)
  - Full pipeline diagram ([PNG](starlight_architecture.png), [PDF](starlight_architecture.pdf), [SVG](starlight_architecture.svg)); shared preprocessing and the three quantification modes (epi2me, isoquant, oarfish).
- [Pipeline overview](pipeline_overview.md)
  - Colleague-facing step-by-step summary (incl. workflow-glue and mode comparison).
- [Implementation status](first_implementations.md)
  - Module inventory, bug-fix log, and open questions.
- [Usage](usage.md)
  - Samplesheet format, example commands per quant mode, barcode extraction parameters.
- [Output](output.md)
  - Published results by pipeline section and quantification mode.
- [Quantification comparison](comparison.md)
  - `compare_quant_methods.py` and `compare_quant_detailed.py` usage.

### Mode-specific design notes

- [IsoQuant implementation plan](isoquant_implementation_plan.md)
- [Oarfish implementation plan](oarfish_implementation_plan.md)
- [QC implementation plan](qc_implementation_plan.md) — read-count funnel, NanoPlot, Seurat (planned)
- [QUIK / barcode implementation](quik_implementations.md)

Repository: [https://github.com/OncoRNALab/panoramaseq_nano](https://github.com/OncoRNALab/panoramaseq_nano)

General nf-core documentation: [https://nf-co.re](https://nf-co.re)
