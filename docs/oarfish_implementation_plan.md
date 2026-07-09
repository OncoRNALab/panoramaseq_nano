# Implementation Plan — Oarfish quantification path

*Status: implemented (v2 — scnanoseq-style FASTQ→transcriptome path; Oarfish 0.9.4; gene aggregation via `oarfish_aggregate_genes.py`).*

---

## Goal

Third quantification mode parallel to **epi2me** and **isoquant**:

| Mode | Param value | Approach |
|------|-------------|----------|
| **epi2me** | `gene_quant_mode = 'epi2me'` | StringTie → txome align → assign_features → create_matrix |
| **isoquant** | `gene_quant_mode = 'isoquant'` | Genome align → UMI dedup → IsoQuant |
| **oarfish** | `gene_quant_mode = 'oarfish'` | R2 FASTQ → txome align → tag → UMI dedup → CB sort → Oarfish → gene aggregation |

Oarfish follows the [nf-core/scnanoseq](https://github.com/nf-core/scnanoseq) transcriptome branch (direct FASTQ alignment, mapped-only filter, dedup after tagging). STARLIGHT uses **Oarfish 0.9.4** (not scnanoseq's 0.6.5) and **CB/UB** tags (not UR).

---

## Architecture

```
QUIK R2 (cDNA FASTQ)
  → gffread(reference GTF) → transcriptome index
  → MINIMAP2_ALIGN (map-ont, -N 100)
  → SAMTOOLS_FILTER_MAPPED (-F 4)
  → TAG_BAM (CB, UB)
  → UMITOOLS_DEDUP
  → SORT_BAM_CB (-t CB)
  → OARFISH (0.9.4)
  → OARFISH_AGGREGATE (gene matrix from reference GTF)
```

Genome alignment is **not** run when `gene_quant_mode = oarfish`.

---

## Modules

| Module | Path | Role |
|--------|------|------|
| `QUANTIFICATION_OARFISH` | `subworkflows/local/quantification_oarfish.nf` | Orchestrates full path |
| `MINIMAP2_ALIGN_TXOME` | nf-core `minimap2/align` | R2 FASTQ → reference transcriptome |
| `SAMTOOLS_FILTER_MAPPED` | `modules/local/samtools_filter_mapped/` | Keep mapped reads only (`-F 4`) |
| `TAG_BAM` | `modules/local/tag_bam/` | Add CB/UB after txome align |
| `SORT_BAM_CB` | `modules/local/sort_bam_cb/` | Collate alignments by cell barcode |
| `OARFISH` | `modules/local/oarfish/` | Single-cell quantification (**v0.9.4**) |
| `OARFISH_AGGREGATE` | `modules/local/oarfish_aggregate/` | Transcript + gene-level 10x MEX output |

Container: `oarfish:0.9.4--h7f5d12c_0`.

---

## Parameters

```groovy
gene_quant_mode                              = 'oarfish'
oarfish_args                                 = '--single-cell --model-coverage --filter-group no-filters'
oarfish_save_transcript_secondary_alignment  = true   // false → minimap2 -N 1
umitools_dedup_method                        = 'directional'
```

Stricter run example: `--oarfish_args '--filter-group nanocount-filters'` and `--oarfish_save_transcript_secondary_alignment false`.

---

## Outputs

Published under `{outdir}/oarfish/{sample_id}/`:

| Output | Description |
|--------|-------------|
| `{sample}_transcript_bc_matrix/` | 10x MEX (transcript × barcode) |
| `{sample}_gene_bc_matrix/` | 10x MEX (gene × barcode; summed from transcripts) |
| `{sample}.aggregation_stats.json` | Feature/barcode counts |
| `{sample}.meta_info.json` | Oarfish run metadata |
| `transcriptome_align/` | Txome-aligned BAM |
| `cb_sorted/` | CB-sorted BAM fed to Oarfish |

---

## Job script example

```bash
nextflow run . -profile singularity \
    --gene_quant_mode oarfish \
    --genome_fasta /path/to/genome.fa \
    --genome_gtf /path/to/genome.gtf \
    --strandedness unstranded \
    ...
```

---

## References

- Oarfish: https://github.com/COMBINE-lab/oarfish
- scnanoseq: https://github.com/nf-core/scnanoseq
- Gene aggregation: https://github.com/COMBINE-lab/oarfish/issues/41
