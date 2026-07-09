# Pipeline Implementation Status — nf-core/starlight

*Last updated: 2026-06-28*

> **Full architecture diagram and mode comparison:** see [architecture.md](architecture.md).

---

## Status summary

| Component | Status |
|-----------|--------|
| Shared preprocessing (Restrander → QUIK → genome align → TAG_BAM) | ✅ implemented |
| epi2me quantification (`gene_quant_mode = 'epi2me'`) | ✅ implemented |
| IsoQuant quantification (`gene_quant_mode = 'isoquant'`) | ✅ implemented |
| Oarfish quantification (`gene_quant_mode = 'oarfish'`) | ✅ implemented |
| MultiQC + version reporting | ✅ implemented |

All three quantification modes share preprocessing through **TAG_BAM**; exactly one branch runs per execution. Quantification outputs publish under **`{outdir}/{gene_quant_mode}/`** so multiple modes can coexist in one outdir.

---

## Quick mode selector

```bash
# epi2me (wf-single-cell style; requires stringtie_enabled)
--gene_quant_mode epi2me --stringtie_enabled true

# IsoQuant on genome BAM (classical long-read quantification)
--gene_quant_mode isoquant \
  --isoquant_args '--data_type nanopore --no_model_construction --complete_genedb ...'

# Oarfish single-cell on transcriptome BAM
--gene_quant_mode oarfish \
  --oarfish_args '--single-cell --filter-group no-filters --model-coverage'
```

---

## Implemented modules / scripts

### Scripts (`bin/`)

| Script | Purpose | Status |
|--------|---------|--------|
| `extract_barcode.py` | Spatial barcode + structured UMI extraction via parasail SW | ✅ |
| `split_reads.py` | Split oriented reads into R1 (barcode) + R2 (cDNA) | ✅ |
| `tag_bam.py` | Add CB/CR/CY/UR/UB/UY to genome-aligned BAM | ✅ |
| `extract_bam_tags.py` | Extract mapq + barcode tag TSVs from tagged genome BAM | ✅ |
| `aggregate_matrix.py` | Aggregate HDF chunks → 10x MEX (epi2me) | ✅ |
| `oarfish_aggregate_genes.py` | Transcript MEX → gene-level MEX via reference GTF (oarfish) | ✅ |
| `workflow-glue` | CLI for assign_features, create_matrix (epi2me) | ✅ |
| `orient_reads.py` | Strand orientation via vsearch | ⚠️ deprecated |

### Local modules (`modules/local/`)

| Module | Status | Notes |
|--------|--------|-------|
| `restrander` | ✅ | Mandatory orientation step |
| `extract_barcode`, `split_reads`, `quik_starsolo` | ✅ | Barcode / UMI front-end |
| `tag_bam`, `filter_primary_bam` | ✅ | BAM tagging and filtering |
| `sort_bam_names`, `extract_bam_tags`, `assign_features`, `create_matrix` | ✅ | epi2me path |
| `isoquant` | ✅ | IsoQuant wrapper |
| `minimap2_align_bam`, `sort_bam_cb`, `oarfish`, `oarfish_aggregate` | ✅ | oarfish path |
| `samtools_faidx`, `samtools_index` | ✅ | Index helpers |

Subworkflows: `quantification_epi2me.nf`, `quantification_isoquant.nf`, `quantification_oarfish.nf`

---

## Bug fixes applied

### 1. RT-adapter orientation naming inverted (`orient_reads.py`) — deprecated module
**Date:** 2026-05-12 · Module removed in favour of Restrander (2026-06).

### 2. QUIK strategy name invalid
**Date:** 2026-05-12 / 2026-05-22 · Default updated to `4_7_mer_gpu_v1`.

### 3. `editdistance` has no `__version__` attribute
**Date:** 2026-05-22 · Fixed in `extract_barcode/main.nf`.

### 4. UMI lost after SPLIT_READS
**Date:** 2026-05-22 · UR/UY embedded in R2 header; minimap2 `-y` copies to BAM.

### 5. Whitelist barcode orientation
**Date:** 2026-05-22 · Whitelist must contain RC of RT-primer barcodes.

### 6. MINIMAP2_ALIGN_TRANSCRIPTOME join mismatch
**Date:** 2026-06-25 · Join on `meta.id` only.

### 7. Name-sorted BAM indexing failure
**Date:** 2026-06-26 · Removed index step from `SORT_BAM_NAMES`.

### 8. `--min_read_coverage` argparse type
**Date:** 2026-06-26 · Changed to `type=float` in assign_features.

### 9. Polars / pyarrow container for assign_features
**Date:** 2026-06-26 · Updated Wave container with pyarrow.

### 10. IsoQuant output mismatch with `--no_model_construction`
**Date:** 2026-06-27 · `corrected_reads.bed.gz` optional; dynamic grouped-output symlinks.

### 11. Oarfish MINIMAP2_ALIGN_BAM BAM-as-FASTQ error
**Date:** 2026-06-28 · Convert dedup BAM via `samtools fastq -T CB,UB,…` before minimap2.

### 12. OARFISH_AGGREGATE channel join arity
**Date:** 2026-06-28 · Removed incorrect 6-argument `.map` after channel joins.

### 13. Oarfish gene aggregation script
**Date:** 2026-06-28 · Execute permission on `oarfish_aggregate_genes.py`; GTF-based transcript→gene mapping (Oarfish outputs plain transcript IDs).

---

## Open questions

1. **Oarfish barcode count** — test sample `barcode05` reports 1 barcode in Oarfish output; investigate CB tag preservation through txome realignment vs IsoQuant grouped outputs (~10k barcodes).
2. **IsoQuant gene DB caching** — pre-built `.db` param not yet exposed; `--complete_genedb` recommended via `--isoquant_args`.
3. **Minimum cDNA length** — default 50 bp; confirm against library insert distribution.
4. **Whitelist RC convention** — consider `--rc_barcode_whitelist` helper parameter.

---

## Further reading

- [architecture.md](architecture.md) — full pipeline diagram and mode comparison
- [isoquant_implementation_plan.md](isoquant_implementation_plan.md)
- [oarfish_implementation_plan.md](oarfish_implementation_plan.md)
- [quik_implementations.md](quik_implementations.md)
