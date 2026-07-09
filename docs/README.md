# nf-core/starlight: Documentation

The nf-core/starlight documentation is split into the following pages:

- [Architecture](architecture.md)
  - Full pipeline diagram, shared preprocessing, and the three quantification modes (epi2me, isoquant, oarfish).
- [Implementation status](first_implementations.md)
  - Module inventory, bug-fix log, and open questions.
- [Usage](usage.md)
  - How to run the pipeline and command-line flags (see also `nextflow_schema.json`).
- [Output](output.md)
  - Published results by pipeline section and quantification mode.

### Mode-specific design notes

- [IsoQuant implementation plan](isoquant_implementation_plan.md)
- [Oarfish implementation plan](oarfish_implementation_plan.md)
- [QC implementation plan](qc_implementation_plan.md) — read-count funnel, NanoPlot, Seurat (planned)
- [QUIK / barcode implementation](quik_implementations.md)

General nf-core documentation: [https://nf-co.re](https://nf-co.re)
