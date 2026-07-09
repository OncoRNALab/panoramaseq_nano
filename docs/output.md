# nf-core/starlight: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to `--outdir`.

Quantification results are published under **`{outdir}/{gene_quant_mode}/`**, where `gene_quant_mode` is `epi2me`, `isoquant`, or `oarfish`. Shared preprocessing outputs stay at the top level so all three methods can reuse the same front-end results when comparing modes with `-resume` on the same outdir.

For a workflow overview see [architecture.md](architecture.md).

---

## Output tree

```
{outdir}/
├── fastqc/                         # shared
├── restrander/{sample}/
├── quik_starsolo/{sample}/
├── minimap2/{sample}/
├── tag_bam/{sample}/
├── pipeline_info/
├── epi2me/                         # epi2me run only
├── isoquant/                       # isoquant run only
└── oarfish/                        # oarfish run only
```

After running all three modes on the same sample (same `--outdir`, different `--gene_quant_mode`, `-resume`), all three method folders coexist.

---

## Shared preprocessing outputs

These directories are created at `{outdir}/` for every run (some are conditional on parameters).

### FastQC

- `fastqc/`
  - `*_fastqc.html`, `*_fastqc.zip`

### Chopper (optional)

- `chopper/<sample>/`
  - `<sample>.fastq.gz` — when `--chopper_enabled true`

### Restrander

- `restrander/<sample>/`
  - `<sample>.restranded.fastq.gz`, `<sample>.stats.json`

### Extract barcode

- `extract_barcode/<sample>/`
  - `<sample>.bc_tags.tsv.gz`

### QUIK / Starsolo

- `quik_starsolo/<sample>/`
  - Filtered R1/R2 FASTQ with corrected barcodes

### Genome alignment

- `minimap2/<sample>/`
  - Genome-aligned, coordinate-sorted BAM + BAI

### Tag BAM

- `tag_bam/<sample>/`
  - `<sample>.tagged.bam` + index — `CB`, `CR`, `CY`, `UB`, `UR`, `UY` tags

---

## Quantification outputs

Paths below are relative to **`{outdir}/{gene_quant_mode}/`**. Only one mode runs per execution.

### epi2me (`{outdir}/epi2me/`)

Requires `--gene_quant_mode epi2me` and `--stringtie_enabled true`.

| Directory | Key files |
|-----------|-----------|
| `filter_primary_bam/<sample>/` | Primary-filtered genome BAM |
| `stringtie/<sample>/` | `<sample>.transcripts.gtf`, abundance tables |
| `transcriptome/<sample>/` | Sample transcriptome FASTA |
| `minimap2_transcriptome/<sample>/` | Transcriptome-aligned BAM |
| `gffcompare/<sample>/` | `<sample>.tmap` |
| `extract_bam_tags/<sample>/` | mapq and barcode tag TSVs |
| `assign_features/<sample>/` | `<sample>.feature_assigns.tsv.zst` |
| `create_matrix/<sample>/` | `{sample}_gene_bc_matrix/`, `{sample}_transcript_bc_matrix/`, `{sample}.matrix_stats.json` |
| `multiqc/` | MultiQC report for this run |

---

### isoquant (`{outdir}/isoquant/`)

Requires `--gene_quant_mode isoquant`.

| Directory | Key files |
|-----------|-----------|
| `filter_primary_bam/<sample>/` | Primary-filtered genome BAM |
| `umitools_dedup/<sample>/` | `<sample>.dedup.bam`, dedup log |
| `<sample>/` | IsoQuant output tables and `isoquant.log` |
| `multiqc/` | MultiQC report for this run |

Typical files under `<sample>/` (prefix = sample id):

| File | Description |
|------|-------------|
| `{sample}.gene_grouped_tag_CB_counts.matrix.mtx` | Gene count matrix (MTX) |
| `{sample}.gene_grouped_tag_CB_counts.barcodes.tsv` | Cell barcodes for gene matrix |
| `{sample}.gene_grouped_tag_CB_counts.features.tsv` | Gene features for gene matrix |
| `{sample}.transcript_grouped_tag_CB_counts.matrix.mtx` | Transcript count matrix (MTX) |
| `{sample}.transcript_grouped_tag_CB_counts.barcodes.tsv` | Cell barcodes for transcript matrix |
| `{sample}.transcript_grouped_tag_CB_counts.features.tsv` | Transcript features for transcript matrix |
| `{sample}.gene_counts.tsv` | Sample-level gene counts (aggregate) |
| `{sample}.read_assignments.tsv.gz` | Per-read assignment table |

Grouped per-barcode outputs use MTX format (`--counts_format mtx`). Read group name (`tag_CB`) is embedded in filenames when `--read_group tag:CB` is used.

---

### oarfish (`{outdir}/oarfish/`)

Requires `--gene_quant_mode oarfish`.

| Directory | Key files |
|-----------|-----------|
| `transcriptome/reference/` | Reference transcriptome FASTA (gffread) |
| `tag_bam/<sample>/` | Tagged transcriptome BAM (CB/UB) |
| `umitools_dedup/<sample>/` | `<sample>.dedup.bam`, dedup log |
| `<sample>/transcriptome_align/` | Txome-aligned BAM from R2 FASTQ |
| `<sample>/cb_sorted/` | `<sample>.cb_sorted.bam` |
| `<sample>/quant/` | Raw Oarfish MEX + `{sample}.meta_info.json` |
| `<sample>/` | `{sample}_transcript_bc_matrix/`, `{sample}_gene_bc_matrix/`, `{sample}.aggregation_stats.json` |
| `multiqc/` | MultiQC report for this run |

Gene-level matrices are built from transcript-level Oarfish output using the reference GTF (`oarfish_aggregate_genes.py`).

---

## Pipeline information (shared)

- `pipeline_info/`
  - `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`
  - `params.json` — reflects the **most recent** run (including `gene_quant_mode`)
  - `samplesheet.valid.csv`
  - `nf_core_starlight_software_mqc_versions.yml`

Per-mode MultiQC reports live under `{outdir}/{gene_quant_mode}/multiqc/` and are not overwritten across mode switches.

See [Nextflow tracing docs](https://www.nextflow.io/docs/latest/tracing.html) for report details.
