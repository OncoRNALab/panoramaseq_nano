# Implementation Plan — Classical quantification (UMI-tools + IsoQuant)

*Status: implemented (v1 — epi2me + isoquant + oarfish branching).*

---

## Goal

Add a **classical** quantification path parallel to the existing **epi2me / wf-single-cell** path.
The user selects one mode per run via a pipeline parameter.

| Mode | Param value | Approach |
|------|-------------|----------|
| **epi2me** (current default) | `gene_quant_mode = 'epi2me'` | StringTie → txome align → gffcompare → `workflow-glue assign_features` → `create_matrix` |
| **classical** | `gene_quant_mode = 'isoquant'` | UMI-tools dedup on genome BAM → IsoQuant per-barcode quantification |
| **oarfish** | `gene_quant_mode = 'oarfish'` | UMI-tools dedup → txome realign → Oarfish sc quantification |

**Oarfish** is implemented — see `pipeline/docs/oarfish_implementation_plan.md`.

---

## Architecture

```
                    … → MINIMAP2_ALIGN → TAG_BAM
                                    │
                    ┌───────────────┴───────────────────┐
                    │     params.gene_quant_mode         │
                    └───────────────┬───────────────────┘
         epi2me                     │              isoquant              oarfish
  (StringTie + create_matrix)       │         (this document)      (see oarfish plan)
                                    │
    FILTER_PRIMARY_BAM              │         FILTER_PRIMARY_BAM    FILTER_PRIMARY_BAM
    STRINGTIE → GFFREAD             │         UMITOOLS_DEDUP        UMITOOLS_DEDUP
    MINIMAP2 txome (R2 FASTQ)        │         ISOQUANT              txome realign → OARFISH
    GFFCOMPARE                      │
    ASSIGN_FEATURES                 │
    CREATE_MATRIX                   │
```

Full diagram: `docs/architecture.md`

---

## Existing assets (already in repo)

| Asset | Path |
|-------|------|
| IsoQuant module | `modules/local/isoquant/main.nf` |
| UMI-tools dedup module | `modules/nf-core/umitools/dedup/main.nf` |
| Genome tagger | `modules/local/tag_bam/` + `bin/tag_bam.py` |
| Primary BAM filter | `modules/local/filter_primary_bam/` |
| FAIDX helper | `modules/local/samtools_faidx/` |

---

## Phase 1 — Parameter & workflow branching

### 1.1 New parameters (`nextflow.config`)

```groovy
gene_quant_mode = 'epi2me'   // 'epi2me' | 'isoquant'

// IsoQuant-specific — see "IsoQuant run modes" for --no_model_construction vs discovery
isoquant_args          = '--data_type nanopore --no_model_construction --complete_genedb'
isoquant_read_group    = 'tag:CB'
isoquant_use_primary   = true    // filter to primary alignments before dedup

// UMI-tools dedup
umitools_dedup_method  = 'directional'
umitools_dedup_args    = '--extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-cell'
```

Add to `nextflow_schema.json` as enum with descriptions.

### 1.2 Validation (`subworkflows/local/utils_nfcore_starlight_pipeline/main.nf`)

- `gene_quant_mode` ∈ `{epi2me, isoquant}`
- Both modes require `genome_fasta` (or `genome_mmi`) + `genome_gtf` when genome alignment runs
- **epi2me**: keep existing `stringtie_enabled` checks
- **isoquant**: auto-disable or warn if `stringtie_enabled` is also true (mutually exclusive downstream)

### 1.3 Refactor `workflows/starlight.nf`

Extract two subworkflows:

```
subworkflows/local/quantification_epi2me.nf
subworkflows/local/quantification_isoquant.nf
```

Main workflow after `TAG_BAM`:

```nextflow
if (params.gene_quant_mode == 'epi2me' && params.stringtie_enabled && params.genome_gtf && params.genome_fasta) {
    QUANTIFICATION_EPI2ME(...)
} else if (params.gene_quant_mode == 'isoquant' && params.genome_gtf && params.genome_fasta) {
    QUANTIFICATION_ISOQUANT(...)
}
```

Move `SAMTOOLS_FAIDX` out of the epi2me-only block if both modes need `.fai`
(or run FAIDX inside each subworkflow — acceptable duplication for one reference).

---

## Phase 2 — Classical path modules

### 2.1 Input BAM

**Source:** `TAG_BAM.out.bam` (coordinate-sorted; tags `CB`, `UR`, `UB`, `CR`, `CY`, `UY`).

**Optional filter:** reuse `FILTER_PRIMARY_BAM` (`-F 0x904`) before dedup when
`isoquant_use_primary = true` (recommended).

### 2.2 UMI-tools dedup

**Module:** `modules/nf-core/umitools/dedup/main.nf`

**Recommended `ext.args` in `conf/modules.config`:**

```groovy
withName: 'UMITOOLS_DEDUP' {
    ext.prefix = { "${meta.id}.dedup" }
    ext.args = {
        [
            '--extract-umi-method=tag',
            '--umi-tag=UB',
            '--cell-tag=CB',
            '--per-cell',
            "--method=${params.umitools_dedup_method}",
            params.umitools_dedup_args ?: '',
        ].join(' ')
    }
    publishDir = [
        path: { "${params.outdir}/${params.gene_quant_mode}/umitools_dedup/${meta.id}" },
        mode: params.publish_dir_mode,
        pattern: '*.{bam,log}'
    ]
}
```

**ONT long-read note:** use coordinate + `--per-cell` dedup before IsoQuant assigns genes
(published ONT sc workflow: UMI-tools dedup → IsoQuant quant). `--per-gene` is only needed
if deduplicating by gene before IsoQuant (requires GN/GX tag on BAM — defer unless benchmarking
shows UMI collisions).

**BAM requirements:**

- Coordinate-sorted + indexed (already true after `MINIMAP2_ALIGN`)
- Reads must have `CB` and `UB` tags (`tag_bam.py` sets `UB = UR` initially; dedup corrects `UB`)

**After dedup:** index deduped BAM with `samtools index` (inline in subworkflow or small wrapper process).

### 2.3 IsoQuant

**Module:** `modules/local/isoquant/main.nf`

**Wire inputs per sample:**

```nextflow
ISOQUANT(
    [ meta, dedup_bam, dedup_bai, genome_fasta, genome_fai, genome_gtf ],
    params.isoquant_read_group   // e.g. "tag:CB"
)
```

**Recommended `ext.args`:**

```groovy
withName: 'ISOQUANT' {
    ext.args = {
        [
            '--barcoded_bam',
            '--barcode_tag CB',
            '--umi_tag UB',
            params.isoquant_args ?: '',
        ].join(' ')
    }
    publishDir = [
        path: { "${params.outdir}/${params.gene_quant_mode}/${meta.id}" },
        mode: params.publish_dir_mode,
        pattern: '*/*/*.{tsv,tsv.gz,gtf,bed.gz,log}'
    ]
}
```

Suggested ONT flags depend on the analysis goal. See
[IsoQuant run modes](#isoquant-run-modes) below for the full decision guide.

**Module note:** local module passes `--read_group $group_category`; supply `tag:CB` as the
`group_category` input value.

## IsoQuant run modes

> **Pipeline defaults (current):** IsoQuant flags are assembled in `conf/modules.config`
> (`sensitive_ont`, `default_ont`, conditional `retain_introns` quant modes). See
> [Default sensitivity (relaxed / permissive defaults)](architecture.md#default-sensitivity-relaxed--permissive-defaults)
> in `docs/architecture.md`. The modes below describe **alternative** configurations via
> `--isoquant_args`.

IsoQuant behaviour is controlled mainly by `--no_model_construction`,
`--complete_genedb`, and whether the gene annotation is supplied as a raw GTF or a
pre-built `.db` file. Choose the mode that matches your biological question.

#### Mode A — Quantification only **(fastest; strict annotation-only)**

Use when you only need expression counts against a **known** annotation and do **not**
need novel genes or novel isoforms. Pass via `--isoquant_args` to override the pipeline
default (`sensitive_ont` model construction).

```bash
--isoquant_args '--data_type nanopore --no_model_construction --complete_genedb --transcript_quantification unique_only --gene_quantification unique_only'
```

| Flag | Role |
|------|------|
| `--no_model_construction` | Skip transcript model discovery; quantify against existing annotation only |
| `--complete_genedb` | Fast GTF→DB conversion for Ensembl/GENCODE (required when passing a GTF) |

This is the **fastest** strict quantification-only configuration (not the current pipeline default).

#### Mode B — Novel genes and novel isoforms **(recommended for discovery)**

Use when you want to discover **novel genes** and **novel isoforms**. In this mode,
**do not** pass `--no_model_construction`.

**Recommended setup:** pre-build the gene database once (see
[Pre-build and reuse the gene database](#2-pre-build-and-reuse-the-gene-database-recommended-for-discovery-runs))
and pass the `.db` file to `--genedb` instead of the raw GTF. This avoids repeating
the costly GTF conversion on every run while allowing full model construction.

```bash
--isoquant_args '--data_type nanopore --complete_genedb --transcript_quantification unique_only --gene_quantification unique_only'
# plus: pass pre-built .db to --genedb (pipeline param isoquant_genedb — not yet implemented)
```

Pre-building the database enables IsoQuant to detect a **larger number of isoforms**
(including truly novel loci) compared to Mode C below.

#### Mode C — Novel isoforms within annotated genes **(alternative discovery mode)**

Use when you want **novel isoform models** but only for genes **already present** in
your GTF — not novel genes at unannotated loci.

```bash
--isoquant_args '--data_type nanopore --complete_genedb --transcript_quantification unique_only --gene_quantification unique_only'
# omit --no_model_construction; pass GTF (or pre-built .db) to --genedb
```

| Aspect | Mode B (pre-built `.db`, no `--no_model_construction`) | Mode C (`--complete_genedb`, no `--no_model_construction`) |
|--------|--------------------------------------------------------|-------------------------------------------------------------|
| Novel genes | Yes | No — limited to genes in the GTF |
| Novel isoforms | Yes; broadest detection | Yes; only for genes already in the GTF |
| GTF→DB cost | Zero if `.db` pre-built | ~10–20 min with `--complete_genedb` (or hours without it) |
| Runtime | Slower than Mode A | Slower than Mode A; fewer isoforms than Mode B |

**Summary**

| Goal | `--no_model_construction` | `--complete_genedb` | Gene annotation input |
|------|---------------------------|---------------------|------------------------|
| Quantification only (fastest) | **Yes** | **Yes** | GTF or pre-built `.db` |
| Novel genes + isoforms | **No** | Yes (at DB build time) | **Pre-built `.db` recommended** |
| Novel isoforms within GTF genes | **No** | **Yes** | GTF or pre-built `.db` |

**Important:** `--complete_genedb` is required for Ensembl/GENCODE GTFs whenever a
GTF is passed to `--genedb`. Without it, gffutils re-infers gene/transcript features
that already exist in the annotation, which can turn GTF→SQLite conversion from
~10–20 minutes into several hours. See
[Performance & future improvements](#performance--future-improvements).

### 2.4 Outputs

| Output | IsoQuant file |
|--------|---------------|
| Per-barcode gene counts | `{prefix}.gene_grouped_counts.tsv` |
| Per-barcode transcript counts | `{prefix}.transcript_grouped_counts.tsv` |
| Read assignments | `{prefix}.read_assignments.tsv.gz` |
| TPM / count tables | `*_grouped_tpm.tsv`, `*_counts.tsv` |

These are **IsoQuant TSV matrices** (features × barcodes), not 10x MEX.
See Phase 4 for optional MEX conversion.

---

## Phase 3 — Containers & resources

| Step | Container |
|------|-----------|
| UMITOOLS_DEDUP | nf-core Wave `umi_tools_future_matplotlib_numpy_pruned` (in module) |
| ISOQUANT | biocontainers `isoquant:3.12.0` (in module; ≥3.12 required for `--barcoded_bam`) |
| FAIDX / index | `pysam_samtools_python` Wave image |

No new Wave image required unless combining dedup + index in one process.

---

## Phase 4 — Optional enhancements

### 4.1 IsoQuant → 10x MEX converter

Script `bin/isoquant_to_mex.py` to convert `*_gene_grouped_counts.tsv` →
`{sample}_gene_bc_matrix/` for Seurat/Scanpy parity with epi2me output.

### 4.2 Oarfish alternative

- Add `modules/local/oarfish/` when module is available
- Same upstream: `TAG_BAM` → optional filter → UMI-tools dedup
- Third enum value `gene_quant_mode = 'oarfish'`

### 4.3 Whitelist filtering

Optional pre-filter: drop reads with `CB` not in `--barcode_whitelist` before dedup.

---

## Phase 5 — Documentation & testing

1. Update `pipeline/docs/first_implementations.md` — dual-path overview
2. Promote this file to `pipeline/docs/isoquant_implementation.md` once implemented
3. Job script examples:

```bash
# epi2me (default)
--gene_quant_mode epi2me --stringtie_enabled true

# isoquant — Mode A: quantification only (fastest)
--gene_quant_mode isoquant \
    --isoquant_args '--data_type nanopore --no_model_construction --complete_genedb --transcript_quantification unique_only --gene_quantification unique_only'

# isoquant — Mode B: novel genes + isoforms (pre-built .db recommended; omit --no_model_construction)
--gene_quant_mode isoquant \
    --isoquant_args '--data_type nanopore --complete_genedb --transcript_quantification unique_only --gene_quantification unique_only'

# isoquant — Mode C: novel isoforms within GTF genes only (omit --no_model_construction)
# same flags as Mode B when passing GTF; see "IsoQuant run modes" section
```

4. **Test plan:**
   - Re-run `barcode05` with `-resume` in `isoquant` mode (upstream cached)
   - Compare cell count, total UMIs, sparsity vs epi2me MEX
   - Validate IsoQuant grouped TSV dimensions vs unique `CB` tags in BAM

---

## Implementation order

| Step | Task | Effort |
|------|------|--------|
| 1 | Add `gene_quant_mode` param + schema + validation | Small |
| 2 | Extract `QUANTIFICATION_EPI2ME` subworkflow (move existing block) | Small |
| 3 | Create `QUANTIFICATION_ISOQUANT` subworkflow | Medium |
| 4 | Wire `FILTER_PRIMARY_BAM` → `UMITOOLS_DEDUP` → index → `ISOQUANT` | Medium |
| 5 | `modules.config` publishDir + ext.args | Small |
| 6 | Test on `barcode05`; tune IsoQuant ONT flags | Medium |
| 7 | Docs update | Small |
| 8 | Optional MEX converter | Small |

---

## Performance & future improvements

Observed on VSC (accelgor, ~327k dedup reads, GRCh38 + chrIS/spike-ins GTF): a 4 h
SLURM walltime was exhausted almost entirely by **GTF→SQLite database conversion**
(~4 h) before read quantification could finish. Read count was not the bottleneck.

Root causes and planned mitigations:

### 1. `--complete_genedb` (implemented in job scripts)

For official annotations (Ensembl, GENCODE), IsoQuant must receive
`--complete_genedb` whenever a GTF is passed to `--genedb`. This sets
`disable_infer_genes` and `disable_infer_transcripts` in gffutils and avoids
redundant feature merging.

- **Status:** passed via `--isoquant_args` in `jobs/FAP_UmiW200_BC02.sh` (Mode A)
- **Expected effect:** GTF→DB drops from hours to ~10–20 minutes
- **Note:** `--complete_genedb` affects GTF conversion speed; it is independent of
  `--no_model_construction`. The latter controls whether IsoQuant discovers new
  transcript models — see [IsoQuant run modes](#isoquant-run-modes).
- **Docs:** https://ablab.github.io/IsoQuant/cmd.html

### 2. Pre-build and reuse the gene database (recommended for discovery runs)

The pipeline currently passes the raw GTF to `--genedb`, so IsoQuant rebuilds the
`.db` file on every run (and once per barcode sample). IsoQuant accepts a
pre-built gffutils database directly.

**One-time build** (run on a login or interactive node with the IsoQuant container):

```bash
GENOME_GTF=/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf
GENEDB_OUT=/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.db

# Inside the IsoQuant container (or after loading a module that provides isoquant.py):
python3 isoquant_lib/gtf2db.py "$GENOME_GTF" "$GENEDB_OUT" --complete_genedb

# Alternatively, copy the .db from a completed IsoQuant work directory:
#   work/<hash>/barcode05/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.db
# IsoQuant logs: "Provide this database next time to avoid excessive conversion"

# Reuse on subsequent runs:
#   --genedb /path/to/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.db
```

**Pipeline change (not yet implemented):**

- Add parameter `isoquant_genedb` (path to pre-built `.db`; overrides GTF when set)
- Add process `PREPARE_ISOQUANT_GENEDB` that runs once per GTF version and caches
  the `.db` in the project cache directory
- Wire `ISOQUANT` to pass `.db` to `--genedb` instead of `.gtf`

**Expected effect:** zero GTF conversion time on subsequent runs. **Required for
Mode B** (novel gene/isoform discovery) to avoid paying the DB conversion cost on
every run while model construction is enabled.

### 3. `--genedb_output` on fast local scratch

SQLite databases cannot be created reliably on some shared/NFS mounts. IsoQuant
supports `--genedb_output` to write the annotation DB to a local path:

```
--genedb_output /scratch/gent/vo/000/gvo00027/${USER}/isoquant_genedb_cache
```

**Pipeline change (not yet implemented):** add `isoquant_genedb_output` param and
pass it through `ext.args` in `conf/modules.config`.

### 4. SLURM walltime vs Nextflow process labels

The pipeline currently runs with the **local executor** inside a single SLURM job.
Changing IsoQuant to `process_high` (`time = 16.h` in `conf/base.config`) does
**not** extend the SLURM walltime — the outer `#PBS`/`#SBATCH` limit still applies
to the entire Nextflow session.

**Recommendations:**

- Use `#SBATCH --time=08:00:00` (or longer) for first runs while the gene DB is
  still being built
- After pre-building the `.db`, 4 h should be sufficient for ~350k reads
- Longer term: add a `-profile slurm` so heavy steps (IsoQuant) get their own
  SLURM allocations with independent walltimes

### 5. Large custom reference (140 chromosomes/contigs)

The project reference includes GRCh38 primary assembly plus ERCC spike-ins, chrIS,
and alt contigs. IsoQuant iterates all chromosomes present in the reference even
when most have few or no reads.

**Options:**

- Use a stripped reference (primary assembly only) when spike-ins / chrIS are not
  needed for quantification
- Accept the overhead when those contigs are required for spike-in normalisation

### 6. Per-barcode grouping at scale

With `--read_group tag:CB` and `--barcoded_bam`, IsoQuant groups counts by every
unique `CB` tag. Datasets with hundreds of thousands of barcodes (even at low
read depth) can slow the grouping step after alignment assignment.

**Options to evaluate after DB conversion is fixed:**

- Monitor runtime with current settings first
- Consider `--barcode2spot` / `barcode_spot` grouping if barcodes map to fewer spots
- Use `--read_group none` only for bulk/debug runs (not for final scRNA output)

### 7. Alternative quantification paths

| Option | When to use |
|--------|-------------|
| **epi2me** (`gene_quant_mode = 'epi2me'`) | Default; faster at current data scale; produces 10x-style MEX |
| **Oarfish** (`gene_quant_mode = 'oarfish'`, planned) | UMI long-read quant; may be faster than IsoQuant for ONT sc |
| **IsoQuant bulk mode** (`--read_group none`) | Debugging only; skips per-barcode grouping |

### 8. Runtime expectations (after fixes)

For ~327k dedup reads on GRCh38:

| Setup | Approximate IsoQuant time |
|-------|---------------------------|
| Mode A: `--no_model_construction` + `--complete_genedb`, rebuild DB each run | ~30–90 min (DB build ~10–20 min with `--complete_genedb`) |
| Mode A without `--complete_genedb` | 4 h+ (observed: walltime exceeded during DB build) |
| Mode A/B/C with pre-built `.db` | ~15–45 min (quantification); longer if model construction enabled |
| Mode B/C with model construction + large barcode count | depends on CB grouping and discovery load; monitor |

### 9. Monitoring running jobs

Editors (including Cursor) do not auto-refresh SLURM log files. Use:

```bash
# Overall Nextflow progress
tail -f slurm-<JOBID>.out

# Current IsoQuant step detail (replace hash from Nextflow log)
tail -f /scratch/gent/vo/000/gvo00027/${USER}/work/<hash>/*/.command.log
```

---

## Open decisions (confirm before coding)

1. **Default mode** — keep `epi2me` as default? *(recommended: yes)*
2. **Primary-only BAM** — filter with `-F 0x904` before dedup? *(recommended: yes)*
3. **IsoQuant model construction** — use Mode A (`--no_model_construction`) for
   quantification-only runs; use Mode B (pre-built `.db`, no `--no_model_construction`)
   for novel gene/isoform discovery; Mode C for novel isoforms within GTF genes only
   (see [IsoQuant run modes](#isoquant-run-modes))
4. **Output format** — IsoQuant TSV only, or also 10x MEX in v1?
5. **Oarfish** — include in v1 or defer to v2?

---

## References

- IsoQuant CLI: https://ablab.github.io/IsoQuant/cmd.html
- IsoQuant single-cell mode: https://ablab.github.io/IsoQuant/single_cell.html
- UMI-tools dedup: https://umi-tools.readthedocs.io/en/latest/reference/dedup.html
- ONT sc comparison (UMI-tools → IsoQuant): bioRxiv 2025.07.03.662955
- Current epi2me path: `pipeline/docs/first_implementations.md`
- Detailed QUIK notes: `pipeline/docs/quik_implementations.md`
