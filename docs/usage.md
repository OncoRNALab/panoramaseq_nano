# nf-core/starlight: Usage

STARLIGHT (also published as [panoramaseq_nano](https://github.com/OncoRNALab/panoramaseq_nano)) processes **single-end Oxford Nanopore (ONT)** long-read spatial / single-cell cDNA data.

Pipeline parameters are defined in [`nextflow_schema.json`](../nextflow_schema.json). For workflow structure see [architecture.md](architecture.md) and [pipeline_overview.md](pipeline_overview.md).

---

## Samplesheet input

Create a CSV samplesheet and pass it with `--input`. The schema requires at least **`sample`** and **`fastq_1`**; `fastq_2` is optional (leave empty for single-end ONT).

```csv title="samplesheet.csv"
sample,fastq_1
barcode05,/path/to/BC05_500K.fastq.gz
```

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier (no spaces; used as output prefix) |
| `fastq_1` | Path to gzipped ONT FASTQ (`.fastq.gz` or `.fq.gz`) |
| `fastq_2` | Optional; not used for standard ONT single-end runs |

Multiple rows with the same `sample` name are concatenated before processing (same datatype only).

Example: [`assets/schema_input.json`](../assets/schema_input.json)

---

## Required inputs

| Parameter | Description |
|-----------|-------------|
| `--input` | Samplesheet CSV |
| `--outdir` | Output directory |
| `--restrander_config` | Restrander JSON (default: `assets/restrander/PCB109.json`) |
| `--barcode_whitelist` | CSV whitelist for QUIK (first column = barcode sequence) |
| `--genome_fasta` | Reference genome FASTA |
| `--genome_gtf` | Reference annotation GTF |

Optional: `--genome_mmi` (prebuilt minimap2 index), `--n_reads` (subsample for testing).

---

## Profiles

On VSC / HPC with Singularity and GPU (QUIK):

```bash
-profile singularity,gpu
```

Other nf-core profiles (`docker`, `singularity`, `conda`, …) apply as usual. QUIK requires a GPU-capable container runtime (`--nv` is set in `conf/modules.config` for the GPU profile).

---

## Example: epi2me (default)

```bash
nextflow run /path/to/pipeline -profile singularity,gpu \
  --gene_quant_mode epi2me \
  --stringtie_enabled true \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/outdir/epi2me_run \
  --restrander_config /path/to/pipeline/assets/restrander/PCB109.json \
  --barcode_whitelist /path/to/barcodes.csv \
  --barcode_length 36 \
  --genome_fasta /path/to/genome.fa \
  --genome_gtf /path/to/annotation.gtf
```

---

## Example: isoquant

Same preprocessing; quantification uses genome BAM + UMI-tools dedup + IsoQuant 3.13.

```bash
nextflow run /path/to/pipeline -profile singularity,gpu \
  --gene_quant_mode isoquant \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/outdir/isoquant_run \
  --restrander_config /path/to/pipeline/assets/restrander/PCB109.json \
  --barcode_whitelist /path/to/barcodes.csv \
  --genome_fasta /path/to/genome.fa \
  --genome_gtf /path/to/annotation.gtf \
  -resume
```

Optional stricter/faster IsoQuant behaviour via `--isoquant_args` (see [isoquant_implementation_plan.md](isoquant_implementation_plan.md)).

---

## Example: oarfish

scnanoseq-style path: QUIK R2 FASTQ → reference transcriptome align → tag → dedup → Oarfish 0.9.4. **No genome alignment** in this mode.

```bash
nextflow run /path/to/pipeline -profile singularity,gpu \
  --gene_quant_mode oarfish \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/outdir/oarfish_run \
  --restrander_config /path/to/pipeline/assets/restrander/PCB109.json \
  --barcode_whitelist /path/to/barcodes.csv \
  --genome_fasta /path/to/genome.fa \
  --genome_gtf /path/to/annotation.gtf \
  -resume
```

---

## Comparing all three modes on one sample

Run three times with the **same `--input`** and reference paths. You may use separate outdirs (simplest) or the same `--outdir` with different `--gene_quant_mode` and `-resume` so shared preprocessing is reused.

Then run the comparison scripts — see [comparison.md](comparison.md).

---

## Barcode extraction parameters

Spatial barcode and UMI extraction (`EXTRACT_BARCODE`) is controlled by the **Barcode localization options** group in `nextflow_schema.json`. Commonly tuned parameters:

| Parameter | Default | Role |
|-----------|---------|------|
| `barcode_length` | 36 | Spatial barcode length |
| `barcode_window` | 150 | 3′ query window for barcode alignment |
| `umi_window` | 120 | 5′ query window for UMI alignment |
| `barcode_min_qv` | 15 | Minimum Phred quality across all barcode bases |
| `barcode_max_anchor_ed` | 5 | Max edit distance for 3′ RT-adapter anchor |
| `max_adapter_ed` | 3 | Max edit distance for 5′ adapter anchor (UMI) |
| `min_barcode_len` | null | Allow shorter barcodes with G-padding (strict if unset) |
| `min_cdna_len` | 50 | Minimum R2 length after split |

Design details: [quik_implementations.md](quik_implementations.md).

---

## Params file

For repeated runs, use `-params-file params.yaml`:

```yaml
input: "/path/to/samplesheet.csv"
outdir: "/path/to/results"
gene_quant_mode: "isoquant"
restrander_config: "/path/to/pipeline/assets/restrander/PCB109.json"
barcode_whitelist: "/path/to/barcodes.csv"
genome_fasta: "/path/to/genome.fa"
genome_gtf: "/path/to/annotation.gtf"
barcode_length: 36
```

```bash
nextflow run /path/to/pipeline -profile singularity,gpu -params-file params.yaml
```

> [!WARNING]
> Do not pass pipeline parameters via `-c custom.config`; use `-params-file` or CLI flags. Use `-c` only for resource tuning and institutional config.

---

## Resume and reproducibility

```bash
-resume
```

Reuses cached tasks when inputs and parameters match. When comparing quant modes, `-resume` avoids rerunning Restrander, QUIK, and genome alignment if those steps are unchanged.

Pin a git revision when running from the repository:

```bash
nextflow run https://github.com/OncoRNALab/panoramaseq_nano -r main ...
```

---

## Output

Published paths are described in [output.md](output.md). Quantification results live under `{outdir}/{gene_quant_mode}/`.

---

## Further reading

| Document | Contents |
|----------|----------|
| [architecture.md](architecture.md) | Full workflow diagram and parameters |
| [comparison.md](comparison.md) | Post-run quant method comparison |
| [output.md](output.md) | Output file reference |
| General nf-core usage | [https://nf-co.re/docs/usage](https://nf-co.re/docs/usage) |
