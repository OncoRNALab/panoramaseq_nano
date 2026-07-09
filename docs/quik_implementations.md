# QUIK Barcode Correction Implementation Plan

## Library structure (reference)

From 5' to 3' of the oriented read:

```
[Outer BC] - [5' Adapter: ...TCTGTTGGTGCTGATATTGC] - [TT VVVV TT VVVV TT VVVV TT VVVV TTT] - [GGG] - [cDNA] - [oligodT] - [Spatial BC 36bp] - [RC(RT-adapter) 45bp] - [Outer BC RC]
                              └── adapter anchor ──┘   └────── Structured UMI 27bp ─────────┘  TSO
                                   (20bp)                      16 variable V positions                    └─── 3' end alignment target ────────────────────────────────────────────┘
```

**Key sequences:**
- 5' adapter anchor (Fw constant):  `TCTGTTGGTGCTGATATTGC`  (20 bp)
- Structured UMI pattern:           `TT VVVV TT VVVV TT VVVV TT VVVV TTT`  (27 bp; V = A/C/G)
- TSO:                              `GGG`  (mGmGmG normalised)
- RT-adapter (3' anchor, forward):  `CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG`  (45 bp)
- RT-adapter RC (in oriented read): `CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG`  (45 bp)
- Spatial barcode length:           36 bp

---

## Pipeline flow

```
CHOPPER
    │ trimmed.fastq.gz
    ▼
RESTRANDER
    │ restranded.fastq.gz
    ▼
EXTRACT_BARCODE  (two SW alignments per read: 5' window for UMI, 3' window for spatial BC)
    │
    ├──→  bc_tags.tsv.gz
    │       columns: read_id, cdna_start, barcode_start, CR, CY, anchor_ed, UR, UY
    │       - cdna_start   = position where cDNA begins (just after TSO, from 5' alignment)
    │       - barcode_start = position of spatial barcode (from 3' alignment)
    │       - CR / CY      = raw spatial barcode + quality
    │       - UR / UY      = raw 16-base UMI + quality
    │
    ▼
SPLIT_READS  (bin/split_reads.py)
    │  Uses cdna_start and barcode_start from bc_tags per read:
    ├──→  {sample}_R1.fastq.gz  =  read[barcode_start : barcode_start+36]
    │                               36 bp spatial barcode, always at position 0
    └──→  {sample}_R2.fastq.gz  =  read[cdna_start : barcode_start]
    │                               clean cDNA only (no adapters, no UMI, no barcode)
    │     Reads absent from bc_tags are dropped.
    │
    ▼
QUIK_STARSOLO  ←── whitelist (barcodes.csv)
    │  barcode_start=0, barcode_length=36
    │  Input: R1 (barcode) + R2 (cDNA) as synchronized pair
    ├──→  {sample}_R1_filtered.fastq.gz  (corrected barcode embedded in read ID)
    └──→  {sample}_R2_filtered.fastq.gz  (cDNA reads with accepted barcodes)
    │  Read ID format: @<original_id>_calledidx_<idx>_<corrected_barcode>
    │
    ▼
MINIMAP2  on R2_filtered only
    │  minimap2 -ax splice -uf --MD --junc-bed ref_genes.bed genome.mmi R2_filtered.fastq.gz
    │  | samtools sort → {sample}.sorted.bam + .bai
    │
    ▼
TAG_BAM  (bin/tag_bam.py)
    │  Join BAM read_ids with bc_tags + parse QUIK-corrected barcode from read name
    │  BAM tags written:
    │    CB = corrected spatial barcode  (from QUIK read ID suffix)
    │    CR = raw spatial barcode        (from bc_tags)
    │    CY = barcode quality            (from bc_tags)
    │    UR = raw UMI                    (from bc_tags)
    │    UB = UMI (used for dedup)       (= UR; full UMI correction optional later)
    │    UY = UMI quality               (from bc_tags)
    │  Output: {sample}.tagged.bam
    │
    ▼
ASSIGN_FEATURES  (bin/assign_features.py)
    │  Strand-aware intersection of aligned reads with gene annotation (GTF)
    │  Output: {sample}.feature_assigns.tsv.gz  (read_id, gene, transcript, mapq)
    │
    ▼
CREATE_MATRIX  (bin/create_matrix.py)
    │  Join bc_tags + feature_assigns on read_id
    │  Per (CB, gene): cluster UMIs with umi-tools directional + Levenshtein distance
    │    → each cluster = 1 molecule → UB tag assigned to representative UMI
    │  Count unique UB per (CB, gene) → sparse count matrix
    │  Output: {sample}_raw_feature_bc_matrix/  (10x MEX format)
    │    ├── barcodes.tsv.gz
    │    ├── features.tsv.gz
    │    └── matrix.mtx.gz
```

---

## Step-by-step implementation details

### 1. EXTRACT_BARCODE — add 5' UMI alignment (update existing)

**File to update:** `bin/extract_barcode.py`

Two SW alignments are run per read:

#### 3' alignment (spatial barcode) — already implemented

- Query: last 150 bp of oriented read
- Probe: `poly-A(18) + N*36 + RC(RT-adapter)`
  = `AAAAAAAAAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG`
- Extracts: `barcode_start`, `CR`, `CY`, `anchor_ed`

#### 5' alignment (UMI) — new

- Query: **first 80 bp** of oriented read
- Probe (50 bp):
  ```
  TCTGTTGGTGCTGATATTGC + TT + NNNN + TT + NNNN + TT + NNNN + TT + NNNN + TTT + GGG
  = TCTGTTGGTGCTGATATTGCTTNNNNTTNNNNTTNNNNTTNNNNTTTGGG
  ```
  - Positions  0–19: adapter anchor `TCTGTTGGTGCTGATATTGC`
  - Positions 20–21: fixed `TT`
  - Positions 22–25: variable `NNNN` (UMI block 1)
  - Positions 26–27: fixed `TT`
  - Positions 28–31: variable `NNNN` (UMI block 2)
  - Positions 32–33: fixed `TT`
  - Positions 34–37: variable `NNNN` (UMI block 3)
  - Positions 38–39: fixed `TT`
  - Positions 40–43: variable `NNNN` (UMI block 4)
  - Positions 44–46: fixed `TTT`
  - Positions 47–49: `GGG` (TSO anchor)

- **Anchor check**: `qry_aln[:bc_start_col]` compared to `TCTGTTGGTGCTGATATTGC`
  (edit distance ≤ `max_adapter_ed`, default 3)

- **UMI extraction**: walk through the gapped reference alignment; for each column where
  `ref_aln[i] == 'N'`, take `qry_aln[i]` (removing gaps) → concatenate 16 bases = `UR`

- **cdna_start**: `p_aln.end_query + 1` from the 5' alignment = first base of cDNA in the
  oriented read (just after the GGG TSO)

- **Adapter anchor check for 5' alignment** (new filter):
  - `adapter_hit = qry_aln[:adapter_end_col].replace("-", "")`
  - `adapter_ed = editdistance.eval(adapter_hit, "TCTGTTGGTGCTGATATTGC")`
  - Reads with `adapter_ed > max_adapter_ed` still output barcode but `UR`/`UY` = `""`

#### Updated bc_tags.tsv.gz columns

| Column | Source | Description |
|--------|--------|-------------|
| `read_id` | read | Read identifier |
| `cdna_start` | 5' alignment | First base of cDNA in oriented read |
| `barcode_start` | 3' alignment | First base of spatial barcode in oriented read |
| `CR` | 3' alignment | Raw spatial barcode sequence (36 bp) |
| `CY` | 3' alignment | Spatial barcode quality string (36 chars) |
| `anchor_ed` | 3' alignment | Edit distance of RC(RT-adapter) match |
| `UR` | 5' alignment | Raw 16-base UMI (V positions only) |
| `UY` | 5' alignment | UMI quality string (16 chars) |

Reads where only the 3' alignment succeeded get `cdna_start=0`, `UR=""`, `UY=""`.
Reads where only the 5' alignment succeeded are dropped (no barcode = unusable).

---

### 2. SPLIT_READS module (new)

**New files:**
- `bin/split_reads.py`
- `modules/local/split_reads/main.nf`
- `modules/local/split_reads/environment.yml`

**Script logic (`bin/split_reads.py`):**
1. Load `bc_tags.tsv.gz` → dict `{read_id → (cdna_start, barcode_start)}`
2. Iterate over `restranded.fastq.gz`
3. For each read present in the dict:
   - **R1**: `seq[barcode_start : barcode_start + barcode_length]`  (36 bp barcode)
   - **R2**: `seq[cdna_start : barcode_start]`  (clean cDNA, no adapters)
   - Drop read if R2 length < `min_cdna_len` (default 50 bp)
4. Write R1 → `{prefix}_R1.fastq.gz`, R2 → `{prefix}_R2.fastq.gz`

**Note:** Using `cdna_start` instead of `0` for R2 removes the outer barcode, 5' adapter,
UMI, and TSO — giving minimap2 a clean transcript sequence.

**Parameters:**
- `--barcode_length 36`
- `--min_cdna_len 50`

**Module inputs:**
```nextflow
tuple val(meta), path(oriented_reads)
tuple val(meta), path(bc_tags)
```

**Module outputs:**
```nextflow
tuple val(meta), path("*_R1.fastq.gz"), emit: r1
tuple val(meta), path("*_R2.fastq.gz"), emit: r2
```

**Conda environment:** python only (no extra deps)

---

### 3. QUIK_STARSOLO module (update existing)

**Changes to `modules/local/quik_starsolo/main.nf`:**
- Replace `val(barcode_start)` input with `path(r1)` and `path(r2)`
- Remove the `cp input_R1.fastq input_R2.fastq` mirror hack
- Hardcode `barcode_start=0` in the QUIK binary call

**New input block:**
```nextflow
input:
tuple val(meta), path(r1), path(r2)
path  barcode_file
val   barcode_length
```

**Key change in script block:**
```bash
gunzip -c ${r1} > input_R1.fastq
gunzip -c ${r2} > input_R2.fastq

${QUIK_EXEC} \
    barcodes_only.txt \
    input_R1.fastq \
    input_R2.fastq \
    0 \
    ${barcode_length} \
    ...
```

---

### 4. MINIMAP2 modules

Two nf-core modules are already installed at `modules/nf-core/minimap2/`:

#### 4a. MINIMAP2_INDEX (use nf-core module as-is)

**Module:** `modules/nf-core/minimap2/index/main.nf`
**Purpose:** Build a `.mmi` index from a reference FASTA once per genome.
**Conda env:** `minimap2=2.30` (already in `environment.yml`)

This module is used **only when no prebuilt index is supplied**. A new `genome_mmi`
parameter (optional, default `null`) controls this branching in the workflow:

```nextflow
// In workflows/starlight.nf
if (params.genome_mmi) {
    ch_index = Channel.fromPath(params.genome_mmi).map { [[id:'genome'], it] }
} else {
    MINIMAP2_INDEX(Channel.fromPath(params.genome_fasta).map { [[id:'genome'], it] })
    ch_index = MINIMAP2_INDEX.out.index
}
```

`ext.args` in `conf/modules.config` for `MINIMAP2_INDEX`:
```groovy
ext.args = '-I 16G'   // raise per-thread memory limit for large genomes
```

#### 4b. MINIMAP2_ALIGN (nf-core module, used as-is)

**Module:** `modules/nf-core/minimap2/align/main.nf` — **no modifications**.

Following nf-core convention, all runtime flags are injected via `ext.args` closures
in `conf/modules.config`. The `--junc-bed` argument takes a reference file that lives
on the shared HPC filesystem; passing its absolute path via `params.junc_bed` in
`ext.args` is the standard nf-core pattern for reference files — no per-sample staging
is needed on a shared filesystem.

**Call in workflow (`workflows/starlight.nf`):**
```nextflow
MINIMAP2_ALIGN(
    ch_r2_filtered,   // tuple val(meta), path(reads)
    ch_index,         // tuple val(meta2), path(reference)  — .mmi or .fasta
    true,             // bam_format
    "bai",            // bam_index_extension
    false,            // cigar_paf_format
    false             // cigar_bam
)
```

**`ext.args` in `conf/modules.config`:**
```groovy
withName: 'MINIMAP2_ALIGN' {
    ext.args  = [
        '-ax splice',
        '-uf',
        '--MD',
        '--cap-kalloc 100m',
        params.junc_bed ? "--junc-bed ${params.junc_bed}" : ''
    ].join(' ').trim()
    publishDir = [ path: "${params.outdir}/minimap2", mode: 'copy', pattern: '*.{bam,bai}' ]
}
withName: 'MINIMAP2_INDEX' {
    ext.args = '-I 16G'
    publishDir = [ path: "${params.outdir}/minimap2/index", mode: 'copy', pattern: '*.mmi',
                   enabled: params.save_reference ?: false ]
}
```

**New parameters in nextflow.config:**
- `genome_fasta`    — path to reference genome FASTA (required unless `genome_mmi` given)
- `genome_mmi`      — path to prebuilt `.mmi` index (optional; skips `MINIMAP2_INDEX`)
- `genome_gtf`      — path to gene annotation GTF (required; used by `ASSIGN_FEATURES`)
- `junc_bed`        — path to splice junction BED for `--junc-bed` (optional; if absent,
                       minimap2 still works but junction accuracy may be lower)
- `save_reference`  — boolean; if true, publish the generated `.mmi` index (default `false`)

**Generating `junc_bed` outside the pipeline (one-time step):**
```bash
paftools.js gff2bed annotation.gtf > ref_genes.bed
```
The resulting file is then passed as `--junc_bed ref_genes.bed` at pipeline launch time.
This avoids adding a `paftools.js` process inside the pipeline for an infrequently
changing reference file.

---

### 5. TAG_BAM module (new)

**New files:**
- `bin/tag_bam.py`
- `modules/local/tag_bam/main.nf`
- `modules/local/tag_bam/environment.yml`

**Script logic (`bin/tag_bam.py`):**
1. Load `bc_tags.tsv.gz` → dict `{read_id → (CR, CY, UR, UY)}`
2. Open sorted BAM with pysam
3. For each read:
   - Parse `CB` from read name: split on `_calledidx_` → last token
   - Look up `CR`, `CY`, `UR`, `UY` from bc_tags
   - Set BAM tags: `CB`, `CR`, `CY`, `UR`, `UB` (= `UR`), `UY`
4. Write tagged BAM

**Conda environment:** python, pysam, samtools

---

### 6. ASSIGN_FEATURES module (new)

Adapted from `wf-single-cell`'s `assign_features.py`.

**New files:**
- `bin/assign_features.py`
- `modules/local/assign_features/main.nf`
- `modules/local/assign_features/environment.yml`

**Script logic:**
- Input: tagged BAM + GTF
- Strand-aware intersection of read alignments with gene features (pyranges or pybedtools)
- Reads not overlapping any gene get a genomic interval label (`chr_start_end`) so they
  are still counted (same approach as wf-single-cell `create_region_name`)
- Output: `{sample}.feature_assigns.tsv.gz`  (read_id, gene, transcript, strand, mapq)

**Conda environment:** python, pyranges or pybedtools, pandas

---

### 7. CREATE_MATRIX module (new)

Directly adapted from `wf-single-cell`'s `create_matrix.py`.

**New files:**
- `bin/create_matrix.py`
- `modules/local/create_matrix/main.nf`
- `modules/local/create_matrix/environment.yml`

#### UMI deduplication method

wf-single-cell uses **`umi-tools` directional clustering with Levenshtein distance**
(instead of the default Hamming), monkey-patched to handle ONT indel errors in UMIs.

The algorithm per `(CB, gene)` group:
1. Count occurrences of each raw UMI string (`UR`)
2. Build a directed graph: add edge `umi1 → umi2` if:
   - `edit_distance(umi1, umi2) ≤ 2`  (Levenshtein, handles indels)
   - `count(umi1) ≥ 2 × count(umi2) − 1`  (directional: high-count absorbs low-count)
3. Each connected component = one original molecule
4. Number of components = expression count for that `(CB, gene)`
5. Each read gets its UMI replaced by the representative (most abundant) UMI → `UB` tag

The threshold of 2 is generous for ONT's higher error rate.

```python
# Core code (from wf-single-cell, adapted):
from umi_tools import UMIClusterer
from editdistance import eval as edit_distance

# Monkey-patch: replace Hamming with Levenshtein
def get_adj_list_directional_lev(self, umis, counts, threshold=2):
    adj_list = {umi: [] for umi in umis}
    for umi1, umi2 in itertools.combinations(umis, 2):
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= counts[umi2] * 2 - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= counts[umi1] * 2 - 1:
                adj_list[umi2].append(umi1)
    return adj_list

UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev

# Per (CB, gene) group:
clusterer = UMIClusterer(cluster_method="directional")
umi_counts = collections.Counter(umis)
clusters = clusterer(umi_counts, threshold=2)
# len(clusters) = number of unique molecules
```

**Script logic:**
1. Read `bc_tags.tsv.gz` (with `CB`, `UR`) + `feature_assigns.tsv.gz` (with `gene`)
2. Join on `read_id`
3. Filter: keep only rows where `CB != '-'` and `UR != '-'`
4. Per chromosome chunk (memory-efficient), run `cluster_dataframe`:
   - Group by `(CB, gene)`, apply Levenshtein directional clustering → set `UB`
5. Build sparse count matrix: rows = genes, columns = barcodes, values = unique `UB` counts
6. Output: 10x MEX directory

```
{sample}_raw_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

**Parameters:**
- `--umi_length 16`  (length of extracted UMI V-bases; for validation)
- `--skip_umi_clustering`  (flag for testing: skip clustering, count unique UMI strings)
- `--ref_interval 1000`  (bp window for reads with no gene assignment)

**Conda environment:** python, umi-tools, editdistance, pandas, scipy, numpy

---

## Changes to existing files

| File | Change |
|------|--------|
| `bin/extract_barcode.py` | Add 5' SW alignment for UMI; output `cdna_start`, `UR`, `UY` columns |
| `modules/local/extract_barcode/main.nf` | Pass `--adapter_anchor`, `--umi_pattern` params to script |
| `workflows/starlight.nf` | Wire `SPLIT_READS` → `QUIK_STARSOLO` → `MINIMAP2` → `TAG_BAM` → `ASSIGN_FEATURES` → `CREATE_MATRIX` |
| `modules/local/quik_starsolo/main.nf` | Replace `val(barcode_start)` with `path(r1), path(r2)`; hardcode `barcode_start=0` |
| `nextflow.config` | Add `adapter_anchor`, `genome_fasta`, `genome_mmi`, `genome_gtf`, `junc_bed`, `min_cdna_len`, `max_adapter_ed` |
| `nextflow_schema.json` | Add new parameters |
| `conf/modules.config` | Add `publishDir` stanzas for new modules |

---

## Reused nf-core modules (no changes needed)

| Module | Usage |
|--------|-------|
| `modules/nf-core/minimap2/index/main.nf` | Build `.mmi` index when `genome_mmi` param is not provided |
| `modules/nf-core/minimap2/align/main.nf` | Align R2 cDNA reads; flags via `ext.args` in `modules.config` |

---

## New files to create

| File | Purpose |
|------|---------|
| `bin/split_reads.py` | Split oriented reads into R1 (barcode) + R2 (clean cDNA) |
| `bin/tag_bam.py` | Write CB/CR/CY/UR/UB/UY BAM tags |
| `bin/assign_features.py` | Gene/transcript assignment from tagged BAM + GTF |
| `bin/create_matrix.py` | UMI deduplication + MEX count matrix output |
| `modules/local/split_reads/main.nf` | Nextflow process for split_reads.py |
| `modules/local/split_reads/environment.yml` | Conda env (python) |
| `modules/local/tag_bam/main.nf` | Nextflow process for tag_bam.py |
| `modules/local/tag_bam/environment.yml` | Conda env (python, pysam, samtools) |
| `modules/local/assign_features/main.nf` | Nextflow process |
| `modules/local/assign_features/environment.yml` | Conda env (python, pyranges, pandas) |
| `modules/local/create_matrix/main.nf` | Nextflow process |
| `modules/local/create_matrix/environment.yml` | Conda env (python, umi-tools, editdistance, scipy, pandas) |

---

## Open questions before implementation

1. **Genome reference**: What genome/GTF should be used for testing? (e.g. human GRCh38)
2. **Minimum cDNA length**: What is a reasonable minimum R2 length to pass to minimap2?
   (suggest 50 bp, but depends on shortest transcript of interest)
3. **QUIK GPU fallback**: QUIK requires a GPU + Singularity. Should the pipeline fall back
   to CPU edit-distance correction (like wf-single-cell `assign_barcodes.py`) when no GPU
   is available?
