# Quantification method comparison

Post-pipeline scripts to compare **epi2me**, **isoquant**, and **oarfish** on the same sample when each mode has been run with a separate `--outdir` (or the same `--outdir` with `-resume` and different `--gene_quant_mode`).

Scripts live in `bin/`:

| Script | Purpose |
|--------|---------|
| `compare_quant_methods.py` | Detection metrics per method (genes/transcripts > 0 / > 10, totals, medians) |
| `compare_quant_detailed.py` | Mapping funnel (Restrander → QUIK → downstream), gene overlap (Ensembl IDs), count correlations, markdown + JSON + plots |

---

## Prerequisites

1. Three completed runs (or one shared preprocessing outdir plus three quant subdirs).
2. Reference GTF with `gene_id` and `gene_name` (for epi2me symbol normalization in detailed comparison).
3. Python 3 with `matplotlib` (detailed plots only).

Expected layout (example):

```
FAP_outdir/
├── epi2me_scnanoseq/      # --gene_quant_mode epi2me
├── isoquant_scnanoseq/    # --gene_quant_mode isoquant
└── oarfish_scnanoseq/     # --gene_quant_mode oarfish
```

Each outdir should contain shared preprocessing (`restrander/`, `quik_starsolo/`, …) and its mode folder (`epi2me/`, `isoquant/`, or `oarfish/`).

---

## Basic comparison

```bash
python3 pipeline/bin/compare_quant_methods.py \
  --epi2me-outdir   FAP_outdir/epi2me_scnanoseq \
  --isoquant-outdir FAP_outdir/isoquant_scnanoseq \
  --oarfish-outdir  FAP_outdir/oarfish_scnanoseq \
  --sample barcode05 \
  --output FAP_outdir/comparison_scnanoseq/barcode05_quant_comparison_summary
```

Writes `{output}.md` and `{output}.json` with per-method detection tables.

Optional: `--barcode <CB>` to restrict to one cell barcode; `--gtf` for oarfish gene aggregation checks.

---

## Detailed comparison (funnel + overlap + plots)

```bash
python3 pipeline/bin/compare_quant_detailed.py \
  --epi2me-outdir   FAP_outdir/epi2me_scnanoseq \
  --isoquant-outdir FAP_outdir/isoquant_scnanoseq \
  --oarfish-outdir  FAP_outdir/oarfish_scnanoseq \
  --sample barcode05 \
  --gtf /path/to/Homo_sapiens.GRCh38.109.gtf \
  --output FAP_outdir/comparison_scnanoseq/barcode05_quant_comparison_detailed
```

### Outputs

| File | Content |
|------|---------|
| `{output}.md` | Human-readable report with embedded plot links |
| `{output}.json` | Full structured report |
| `plots/` (default: sibling of output prefix) | Scatter plots, detection bars, gene overlap, mapping funnel |

### Mapping funnel

The funnel plot includes:

- **Restrander** — total oriented reads (upstream)
- **QUIK_STARSOLO (R2 filtered)** — input to downstream alignment/quantification
- Downstream stages labelled as **% of QUIK** (and % of restrander in markdown tables)

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--plots-dir` | `{output_dir}/plots` | PNG output directory |
| `--no-plots` | off | Skip matplotlib figures |
| `--barcode` | all barcodes | Single cell barcode filter |

---

## Interpretation notes

- **Count totals are not directly comparable** across methods without accounting for default sensitivity (Oarfish `no-filters`, IsoQuant `sensitive_ont`, epi2me assignment rules). See [architecture.md](architecture.md#default-sensitivity-relaxed--permissive-defaults).
- **Gene overlap** normalizes epi2me symbols to Ensembl `gene_id` via the reference GTF; unmapped symbols are excluded.
- **Correlations** (Spearman/Pearson) are computed on shared detected genes (count > 0) and are indicative only — count types differ (UMI clustering vs IsoQuant assignment vs Oarfish EM).

See also [output.md](output.md) for published file paths per mode.
