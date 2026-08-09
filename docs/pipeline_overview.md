# Pipeline overview

*summary of pipeline structure and purpose.*

Nextflow pipeline for **Oxford Nanopore (ONT) spatial / single-cell long-read cDNA** data. It takes raw FASTQ reads, extracts cell barcodes and UMIs, aligns cDNA to a reference genome, and produces a **gene/transcript count matrix** per cell.

You choose **one** quantification method per run via `--gene_quant_mode`: **epi2me**, **isoquant**, or **oarfish**.

![STARLIGHT pipeline architecture](starlight_architecture.svg)

See also [architecture.md](architecture.md) for the detailed mermaid diagram and parameter tables.

---

## What goes in and what comes out

**Input:** A samplesheet of single-end ONT FASTQ files (`.fastq.gz`), plus reference genome FASTA and GTF annotation.

**Output:** Per-sample count matrices and QC reports under `--outdir`. Shared preprocessing (alignment, barcode tagging) is written at the top level; quantification results go under `{outdir}/{gene_quant_mode}/`.

---

## Pipeline flow (high level)

```
Raw FASTQ
  → Quality control & optional trimming
  → Read orientation (Restrander)
  → Barcode/UMI extraction & correction
  → Genome alignment
  → BAM tagging (cell + UMI metadata)
  → [ONE OF] epi2me / isoquant / oarfish quantification
  → MultiQC report
```

Everything up to **BAM tagging** is shared across all three quantification modes.

---

## Part 1 — Shared preprocessing

These steps run for every pipeline execution, regardless of quantification mode.

| Step | What it does | Why it's needed |
|------|--------------|-----------------|
| **FASTQC** | Quality control on raw reads | Baseline QC before any processing |
| **SEQTK_SAMPLE** *(optional)* | Subsamples reads to `--n_reads` | Faster testing on a subset of data |
| **CHOPPER** *(optional)* | Quality trimming and length filtering | Removes low-quality bases/reads before downstream steps (`--chopper_enabled true`) |
| **RESTRANDER** | Orients reads to a consistent 5′→3′ direction | ONT reads can arrive in either orientation; all later steps assume a fixed layout (barcode → UMI → cDNA) |
| **EXTRACT_BARCODE** | Finds spatial cell barcode and structured UMI in each read using parasail alignment | Locates where the barcode and UMI sit in the read sequence; records coordinates and quality scores |
| **SPLIT_READS** | Splits each read into **R1** (36 bp barcode) and **R2** (cDNA) | Separates barcode metadata from the cDNA portion that will be aligned and quantified |
| **QUIK_STARSOLO** | GPU-accelerated barcode correction against a whitelist | Corrects sequencing errors in cell barcodes; filters out reads with invalid/unknown barcodes; keeps R1/R2 synchronized |
| **MINIMAP2_INDEX** | Builds a genome index (`.mmi`) | Required for alignment; skipped if a prebuilt index is provided (`--genome_mmi`) |
| **MINIMAP2_ALIGN** | Aligns **R2 cDNA** to the reference genome (`-ax splice`) | Maps each read to a genomic location; UMI information is carried through from the FASTQ header |
| **TAG_BAM** | Adds cell and UMI tags to the aligned BAM | Writes standard BAM tags: **CB** (corrected cell barcode), **UB** (corrected UMI), plus raw/quality tags (CR, CY, UR, UY). All quantification modes depend on these tags |

**Key output of preprocessing:** `{sample}.tagged.bam` — genome-aligned reads with cell barcode and UMI metadata attached.

---

## Part 2 — Quantification (choose one mode)

Only **one** branch runs per execution, selected by `--gene_quant_mode`.

### Option A — **epi2me** (default)

Adapted from Oxford Nanopore's **wf-single-cell** workflow. Builds a **sample-specific transcriptome** and counts UMIs using **workflow-glue** (see below).

| Step | What it does | Why it's needed |
|------|--------------|-----------------|
| **FILTER_PRIMARY_BAM** | Keeps primary alignments only | Removes secondary/supplementary/unmapped reads before assembly |
| **STRINGTIE** | Assembles a per-sample transcriptome from genome alignments | Discovers which isoforms are present in this sample |
| **GFFREAD** | Converts StringTie GTF to transcriptome FASTA | Creates a reference for transcript-level alignment |
| **MINIMAP2 (transcriptome)** | Realigns QUIK-filtered R2 FASTQ to the sample transcriptome | Maps reads directly to isoforms rather than the genome |
| **GFFCOMPARE** | Compares StringTie transcripts to the reference annotation | Produces a mapping file (`.tmap`) linking sample transcripts to known genes |
| **EXTRACT_BAM_TAGS** | Pulls CB/UB tags from the genome BAM | Provides cell/UMI metadata for feature assignment |
| **ASSIGN_FEATURES** *(workflow-glue)* | Assigns each read to gene and transcript features | Links aligned reads to genomic features using alignment quality and coverage thresholds |
| **CREATE_MATRIX** *(workflow-glue)* | Directional UMI clustering → 10x-format count matrix | Deduplicates UMIs per cell and produces gene + transcript count matrices (MEX format) |
| **aggregate_matrix.py** | Converts HDF chunks from create_matrix to 10x MEX | Final packaging of sparse count matrices |

**Best for:** Sample-specific isoform discovery combined with UMI-based counting (Oxford Nanopore wf-single-cell style).

---

### What is workflow-glue? (epi2me only)

**workflow-glue** (also called **wf-glue**) is a Python command-line toolkit **vendored from Oxford Nanopore's [wf-single-cell](https://github.com/nanoporetech/wf-single-cell) pipeline**. In our pipeline it lives under `pipeline/bin/workflow_glue/` and is invoked via the `workflow-glue` CLI wrapper.

It is **not** used by the isoquant or oarfish branches — only epi2me depends on it for the final gene/transcript counting steps.

#### Two commands used in our custom pipeline

| Command | Nextflow module | Role |
|---------|-----------------|------|
| `workflow-glue assign_features` | `ASSIGN_FEATURES` | Feature assignment |
| `workflow-glue create_matrix` | `CREATE_MATRIX` | UMI deduplication + count matrix |

#### `assign_features` — link reads to genes and transcripts

**Inputs:**
- Transcriptome-aligned BAM (reads mapped to the StringTie-assembled transcriptome)
- gffcompare `.tmap` file (maps sample transcript IDs to reference genes)
- Reference GTF (gene names and transcript IDs)
- Tag file from `EXTRACT_BAM_TAGS` (genomic MAPQ per read)

**What it does:**
- Parses transcriptome alignments in chunks
- Applies quality filters: minimum MAPQ (`--min_mapq`, default 30), minimum transcript coverage (`--min_tr_coverage`, default 0.4), minimum read coverage (`--min_read_coverage`, default 0.4)
- Resolves multimapping reads using alignment scores and coverage rules
- Outputs `{sample}.feature_assigns.tsv.zst` — one row per read with assigned gene and transcript

**Why it's needed:** Genome alignment alone does not tell you which isoform a read came from. Realignment to the sample transcriptome plus gffcompare linking allows precise gene/transcript assignment before counting.

#### `create_matrix` — UMI deduplication and count matrix construction

**Inputs:**
- Barcode/UMI tag file from `EXTRACT_BAM_TAGS` (CB, UB, alignment coordinates)
- Feature assignments from `assign_features`

**What it does:**
1. Joins cell barcode + UMI tags with gene/transcript assignments per read
2. Runs **directional UMI clustering** (via UMI-tools `UMIClusterer`, with Levenshtein-distance patching for ONT UMIs) to collapse sequencing errors within each (cell, gene) group
3. Writes intermediate HDF matrix chunks
4. `aggregate_matrix.py` (companion script) merges chunks into final **10x Genomics MEX format** matrices:
   - `{sample}_gene_bc_matrix/` (genes × cells)
   - `{sample}_transcript_bc_matrix/` (transcripts × cells)
5. Writes `{sample}.matrix_stats.json` (barcode counts, feature counts, etc.)

**Why it's needed:** Raw read counts overestimate molecular abundance because of PCR duplication. UMI clustering collapses reads sharing the same biological molecule (allowing for sequencing error), producing integer UMI counts suitable for downstream single-cell analysis (Seurat, Scanpy, etc.).

#### UMI deduplication across all three methods

All three modes use the **directional** UMI deduplication strategy by default, but the **tool, timing, and grouping differ**:

| | **epi2me** | **isoquant** | **oarfish** |
|---|------------|--------------|-------------|
| **Tool** | workflow-glue `create_matrix` (UMI-tools `UMIClusterer` library) | UMI-tools **`dedup` CLI** (`UMITOOLS_DEDUP` module) | UMI-tools **`dedup` CLI** (same module) |
| **Method** | `directional` | `directional` (`--umitools_dedup_method`, default) | `directional` (same param) |
| **When** | After gene/transcript assignment | On tagged genome BAM, before IsoQuant | On tagged **transcriptome** BAM, before Oarfish |
| **Grouping** | Per **(cell barcode, gene)** | Per **cell barcode + alignment position** (`--per-cell` on BAM) | Per cell barcode + alignment position (`--per-cell`) |
| **ONT tweak** | Levenshtein distance (patched in workflow-glue) | Standard UMI-tools `dedup` behaviour | Standard UMI-tools `dedup` behaviour |

Isoquant and oarfish use the same UMI-tools module and `--umitools_dedup_method`, but on different BAMs (genome vs transcriptome). Epi2me skips BAM-level dedup entirely and deduplicates later during matrix construction.

Change the method for isoquant and oarfish together with `--umitools_dedup_method` (e.g. `adjacency`, `directional-adjacency`).

#### How workflow-glue differs from isoquant/oarfish UMI handling

| Aspect | workflow-glue (epi2me) | isoquant / oarfish |
|--------|------------------------|---------------------|
| When dedup runs | After feature assignment, in `create_matrix` | Before quantification, via **UMI-tools dedup** on BAM |
| Dedup scope | Per (cell barcode, gene) | Per cell barcode + genomic position (`--per-cell`) |
| Count unit | Deduplicated UMI clusters (integer) | Read assignments (IsoQuant) or EM estimates (Oarfish) |
| Alignment used | Sample StringTie transcriptome | Genome (IsoQuant) or reference transcriptome (Oarfish) |

---

### Option B — **isoquant**

Quantifies directly from the **genome-aligned BAM** using IsoQuant, with per-cell grouping via the CB tag.

| Step | What it does | Why it's needed |
|------|--------------|-----------------|
| **FILTER_PRIMARY_BAM** | Keeps primary alignments only | Reduces noise from multimapping/secondary alignments before deduplication |
| **UMI-tools dedup** | Collapses PCR duplicates on the tagged genome BAM before IsoQuant | Same UMI-tools module/flags as oarfish; applied on genome BAM |
| **SAMTOOLS_INDEX** | Indexes the deduplicated BAM | Required for IsoQuant input |
| **ISOQUANT** | Long-read quantification on genome BAM with reference annotation | Assigns reads to genes/transcripts using ONT-aware splice correction and model construction; outputs per-cell grouped count tables |

**Default settings** are tuned for sensitivity (`sensitive_ont` model construction, permissive gene quantification). Transcript counts use `unique_only` by default (ambiguous reads excluded at isoform level).

**Best for:** Reference-annotation-based quantification with IsoQuant's ONT-specific read assignment, gene + transcript tables, and optional isoform discovery.

---

### Option C — **oarfish**

scnanoseq-style path: aligns **cDNA FASTQ (R2)** directly to the reference transcriptome, then quantifies with Oarfish 0.9.4. **No genome alignment** is performed in this mode.

| Step | What it does | Why it's needed |
|------|--------------|-----------------|
| **GFFREAD** | Builds transcriptome FASTA from reference GTF | Oarfish requires transcriptome-aligned reads |
| **MINIMAP2_ALIGN_TXOME** | Aligns R2 FASTQ to transcriptome (`map-ont`, `-N 100`) | scnanoseq transcript branch equivalent |
| **SAMTOOLS_FILTER_MAPPED** | Keeps mapped reads only (`-F 4`) | Drop unmapped before tagging |
| **TAG_BAM** | Adds CB/UB tags after txome align | Barcodes applied post-alignment (scnanoseq order) |
| **UMI-tools dedup** | Collapses PCR duplicates on tagged txome BAM | Oarfish counts read records, not raw UMIs |
| **SORT_BAM_CB** | Sorts BAM by cell barcode (`samtools sort -t CB`) | Oarfish single-cell mode requires adjacent CB records |
| **OARFISH** | Probabilistic transcript quantification per cell (v0.9.4) | EM allocation of multimapping reads |
| **OARFISH_AGGREGATE** | Sums transcript counts to gene level | Maps transcripts → genes via reference GTF |

**Best for:** Transcriptome-based single-cell quantification aligned with published scnanoseq benchmarks; typically yields higher counts than isoquant.

---

## Part 3 — Reporting

| Step | What it does | Why it's needed |
|------|--------------|-----------------|
| **MULTIQC** | Aggregates QC metrics into one HTML report | Single dashboard for FastQC, tool versions, and workflow parameters |
| **pipeline_info/** | Nextflow trace, timeline, params snapshot | Reproducibility and debugging |

---

## How the three modes compare

| | **epi2me** | **isoquant** | **oarfish** |
|---|------------|--------------|-------------|
| **Aligns to** | Sample transcriptome (StringTie) | Genome (native) | Reference transcriptome (gffread) |
| **UMI dedup** | workflow-glue `create_matrix`: **directional** clustering, per (CB, gene), after assignment | UMI-tools **`dedup`**: **directional**, `--per-cell`, on BAM before quantification | UMI-tools **`dedup`** on txome BAM after tagging |
| **Multimapping** | Assignment rules in `assign_features` | Discrete assignment classes; `unique_only` at transcript level | EM probabilistic allocation |
| **Isoform discovery** | StringTie assembly | IsoQuant `sensitive_ont` | Reference annotation only |
| **Output format** | 10x MEX (gene + transcript) | IsoQuant TSV / MTX (per cell barcode) | 10x MEX (transcript + gene aggregated) |
| **Sensitivity** | Moderate (sample txome) | Moderately relaxed at gene level; **stricter at transcript level** (`unique_only`) | **Most permissive** (`no-filters` + EM) |

**Important:** Total count numbers are **not directly comparable** across methods — they use different feature universes, count definitions (UMI clusters vs read assignments vs EM estimates), and filtering strictness.

### Critical difference: Oarfish vs IsoQuant (quantification philosophy)

These two methods share the same UMI dedup input but diverge sharply in how reads become transcript counts:

| | **Oarfish** | **IsoQuant** |
|---|-------------|--------------|
| **Default filtering** | `--filter-group no-filters` — disables alignment pre-filters; **does not discard** multimapping/ambiguous alignments before quantification | **Assignment-class filtering** — reads labelled ambiguous, inconsistent, intergenic, etc. are excluded from transcript counts under `unique_only` |
| **Multimapping** | **EM (Expectation–Maximization)** — probabilistically distributes each read across all matching transcripts as **fractional counts** | **Discrete assignment** — each read gets one class; only `unique` (and gene-level `unique_inconsistent`) contribute to counts |
| **Transcript-level default** | Permissive — partial counts spread across all candidate isoforms | Moderately permissive at gene level; **conservative at transcript level** (`transcript_quantification unique_only`) |
| **Philosophy** | **Maximum sensitivity** — catches more transcripts, including lower-confidence assignments | **Maximum specificity** — prioritises high-confidence transcript assignments |

**Observed on `barcode05` (transcripts detected, count > 0):**

- **Oarfish: 22,715** — permissive EM allocates fractional signal across multi-mapped transcripts against the full reference transcriptome (~253k features); txome realignment retains secondaries (`minimap2 -N 100`).
- **IsoQuant: 10,789** — only transcripts with unique or high-confidence assignments enter the reference count matrix; ~172k of 320k assigned reads are excluded at transcript level by assignment filtering.

The higher Oarfish counts are **not necessarily “better”** — they reflect a different quantification philosophy. Choose **Oarfish** when downstream analysis prioritises **sensitivity** (detect as many transcripts as possible, accepting lower-confidence signal). Choose **IsoQuant** when you prioritise **specificity** (count only high-confidence transcript assignments). Also consider **epi2me** for sample-specific isoform discovery with integer UMI counts on a StringTie transcriptome (~14k features).

**Additional factor (Oarfish vs epi2me):** reference transcriptome size (~253k vs ~14k StringTie features) further increases Oarfish’s detected transcript count relative to epi2me — but for **Oarfish vs IsoQuant**, both use the reference annotation; the gap is driven mainly by filtering and EM vs assignment-class logic, not catalogue size alone.

---

## Output directory layout

```
{outdir}/
├── fastqc/                  # Raw read QC
├── restrander/              # Oriented reads
├── quik_starsolo/           # Barcode-corrected R1/R2
├── minimap2/                # Genome-aligned BAM
├── tag_bam/                 # Tagged genome BAM (shared endpoint)
├── pipeline_info/           # Run metadata
│
├── epi2me/                  # ← if --gene_quant_mode epi2me
├── isoquant/                # ← if --gene_quant_mode isoquant
└── oarfish/                 # ← if --gene_quant_mode oarfish
```

Because outputs are scoped by mode, you can run all three methods against the same `--outdir` with `-resume` — shared preprocessing is computed once and reused.

---

## Key parameters to know

| Parameter | Purpose |
|-----------|---------|
| `--gene_quant_mode` | `epi2me` \| `isoquant` \| `oarfish` |
| `--genome_fasta` / `--genome_gtf` | Reference genome and annotation (required) |
| `--barcode_whitelist` | Valid cell barcodes for QUIK correction |
| `--restrander_config` | Read orientation config (e.g. PCB109) |
| `--chopper_enabled` | Enable quality trimming before Restrander |
| `--n_reads` | Subsample for quick test runs |
| `--retain_introns` | IsoQuant only: include intronic/ambiguous reads in counts |
| `--umitools_dedup_method` | IsoQuant + Oarfish: UMI-tools dedup method (default `directional`) |

---

## Quantification comparison — `barcode05` test sample

Three-method comparison on shared preprocessing (`BC05_500K.fastq.gz`, one detected cell barcode). Full tables: `FAP_outdir/new_pipeline/comparison/barcode05_quant_comparison_{wide,long,summary.json,summary.md}`.

### Results (2026-06-28, isoquant re-run with published-style defaults)

| Metric | epi2me | isoquant | oarfish |
|--------|--------|----------|---------|
| Genes detected (>0) | 9,921 | 10,543 | **13,716** |
| Genes detected (>10) | 3,074 | 2,018 | **3,215** |
| Transcripts detected (>0) | 13,579 | 10,789 | **22,715** |
| Transcripts detected (>10) | 2,715 | 1,679 | **3,223** |
| Total gene counts | 287,163 | 204,774 | **321,875** |
| Total transcript counts | 254,933 | 170,341 | **321,875** |

### IsoQuant re-run (slurm-15709759)

Pipeline completed in **~42 min** (ISOQUANT ~40 min). Confirmed command:

`--complete_genedb --stranded none --gene_quantification unique_inconsistent --transcript_quantification unique_only --splice_correction_strategy default_ont --model_construction_strategy sensitive_ont --counts_format mtx`

Compared with the previous strict `--no_model_construction` + `unique_only` gene run:

- **Reference grouped count totals unchanged** (204,774 genes / 170,341 transcripts) — same UMI dedup input (cached) and same assignment buckets; `unique_only` still filters transcript-level counts.
- **New: isoform discovery** — IsoQuant reported **522 novel transcript models** (171 novel-in-catalog, 351 novel-not-in-catalog) plus `transcript_models.gtf` and `extended_annotation.gtf`. These are not reflected in the reference grouped comparison table above.
- **Read funnel:** 349,281 primary → 326,837 dedup → 319,814 assigned (147,309 unique).

### Summary — which method to choose?

| Priority | Recommended method |
|----------|-------------------|
| **Sensitivity** — detect as many transcripts as possible (accept lower confidence) | **Oarfish** |
| **Specificity** — count only high-confidence transcript assignments | **IsoQuant** |
| Sample-specific isoform discovery + integer UMI matrix | **epi2me** |
| Reference annotation + novel IsoQuant GTF outputs | **IsoQuant** (inspect `transcript_models.gtf`, not reference matrix totals alone) |

Cross-method totals are **not directly comparable** (UMI clusters vs read assignments vs EM fractional estimates; different feature universes).

---

## Further reading

| Document | Contents |
|----------|----------|
| [architecture.md](architecture.md) | Detailed technical architecture, module inventory, default sensitivity notes |
| [../FAP_outdir/new_pipeline/comparison/barcode05_quant_comparison_summary.md](../FAP_outdir/new_pipeline/comparison/barcode05_quant_comparison_summary.md) | Three-method quant comparison (`barcode05`) |
| [output.md](output.md) | Published output file reference |
| [usage.md](usage.md) | How to run the pipeline |

Workflow entry: `workflows/starlight.nf`
