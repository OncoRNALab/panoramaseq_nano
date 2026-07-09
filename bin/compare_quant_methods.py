#!/usr/bin/env python3
"""Compare epi2me, isoquant, and oarfish quantification results for one sample.

Metrics (per method, per cell barcode):
  - genes / transcripts detected with count > 0 and count > 10
  - total gene and transcript UMI/count sums
  - number of barcodes in the matrix
  - median count among detected features (count > 0)

Expects the method-scoped output layout: {outdir}/{method}/...
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import statistics
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


METHODS = ("epi2me", "isoquant", "oarfish")


@dataclass
class MethodMetrics:
    method: str
    sample: str
    barcode: str
    n_barcodes: int
    genes_gt_0: int
    genes_gt_10: int
    transcripts_gt_0: int
    transcripts_gt_10: int
    total_gene_counts: float
    total_transcript_counts: float
    median_gene_count_detected: float
    median_transcript_count_detected: float
    source_gene: str
    source_transcript: str
    status: str = "ok"
    note: str = ""


def read_lines(path: Path) -> List[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


def parse_gtf_transcript_gene(gtf_path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = parts[8]
            tx = re.search(r'transcript_id "([^"]+)"', attrs)
            gene = re.search(r'gene_id "([^"]+)"', attrs)
            if tx and gene:
                mapping[tx.group(1)] = gene.group(1)
    return mapping


def read_mtx(path: Path) -> Tuple[int, int, List[Tuple[int, int, float]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                nrows, ncols, _ = map(int, line.split())
                break
        entries = []
        for line in handle:
            row, col, val = line.split()
            entries.append((int(row) - 1, int(col) - 1, float(val)))
    return nrows, ncols, entries


def normalize_barcode(name: str) -> str:
    return name.removesuffix("-1").strip()


def pick_barcode_index(barcodes: List[str], barcode: Optional[str]) -> int:
    if len(barcodes) == 1:
        return 0
    if barcode is None:
        raise ValueError(
            f"Matrix has {len(barcodes)} barcodes; pass --barcode. "
            f"Available: {', '.join(barcodes[:5])}{'...' if len(barcodes) > 5 else ''}"
        )
    target = normalize_barcode(barcode)
    for idx, bc in enumerate(barcodes):
        if normalize_barcode(bc) == target:
            return idx
    raise ValueError(f"Barcode {barcode!r} not found in matrix barcodes")


def counts_for_column(
    nfeatures: int,
    entries: Iterable[Tuple[int, int, float]],
    col_idx: int,
) -> List[float]:
    counts = [0.0] * nfeatures
    for row, col, val in entries:
        if col == col_idx:
            counts[row] = val
    return counts


def detection_metrics(counts: List[float]) -> Dict[str, float]:
    detected = [c for c in counts if c > 0]
    return {
        "detected_gt_0": len(detected),
        "detected_gt_10": sum(1 for c in counts if c > 10),
        "total_counts": sum(counts),
        "median_detected": statistics.median(detected) if detected else 0.0,
    }


def load_mex_counts(matrix_dir: Path, barcode: Optional[str]) -> Tuple[List[float], List[str], str]:
    features_path = matrix_dir / "features.tsv.gz"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    for path in (features_path, barcodes_path, matrix_path):
        if not path.exists():
            raise FileNotFoundError(path)

    feature_ids = []
    for line in read_lines(features_path):
        parts = line.split("\t")
        feature_ids.append(parts[1] if len(parts) >= 2 else parts[0])

    barcodes = read_lines(barcodes_path)
    nrows, ncols, entries = read_mtx(matrix_path)
    col_idx = pick_barcode_index(barcodes, barcode)
    counts = counts_for_column(nrows, entries, col_idx)
    if len(counts) != len(feature_ids):
        raise ValueError(f"Feature/count length mismatch in {matrix_dir}")
    return counts, feature_ids, str(matrix_dir)


def aggregate_transcript_to_gene(
    tx_counts: List[float],
    tx_ids: List[str],
    tx_to_gene: Dict[str, str],
) -> List[float]:
    gene_counts: Dict[str, float] = defaultdict(float)
    for tx_id, val in zip(tx_ids, tx_counts):
        gene_id = tx_to_gene.get(tx_id, tx_id)
        gene_counts[gene_id] += val
    return list(gene_counts.values())


def load_isoquant_grouped(path: Path, barcode: Optional[str]) -> Tuple[List[float], str]:
    with open(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) == 2 and header[0] == "feature_id":
            counts = [float(row[1]) for row in reader if len(row) >= 2]
            return counts, str(path)
        if len(header) < 2:
            raise ValueError(f"Expected grouped or feature_id/count TSV in {path}")
        barcode_cols = header[1:]
        if len(barcode_cols) == 1:
            col_idx = 0
        else:
            col_idx = pick_barcode_index(barcode_cols, barcode)
        counts = []
        for row in reader:
            counts.append(float(row[col_idx + 1]))
    return counts, str(path)


def resolve_isoquant_mtx_prefix(outdir: Path, sample: str, kind: str) -> Optional[Path]:
    base = outdir / "isoquant" / sample
    mtx_path = find_first(
        base,
        [
            f"**/{sample}.{kind}_grouped*counts.matrix.mtx",
            f"**/{sample}.{kind}_grouped*counts.matrix.mtx.gz",
            f"**/{kind}_grouped*counts.matrix.mtx",
        ],
    )
    if mtx_path is None:
        return None
    for suffix in (".matrix.mtx.gz", ".matrix.mtx"):
        if mtx_path.name.endswith(suffix):
            return mtx_path.parent / mtx_path.name[: -len(suffix)]
    return None


def load_isoquant_mtx_counts(prefix: Path, barcode: Optional[str]) -> Tuple[List[float], List[str], List[str], str]:
    matrix_path = None
    for suffix in (".matrix.mtx", ".matrix.mtx.gz"):
        candidate = Path(f"{prefix}{suffix}")
        if candidate.exists():
            matrix_path = candidate
            break
    if matrix_path is None:
        raise FileNotFoundError(f"Missing IsoQuant matrix file for prefix {prefix}")

    barcodes_path = None
    for suffix in (".barcodes.tsv", ".barcodes.tsv.gz"):
        candidate = Path(f"{prefix}{suffix}")
        if candidate.exists():
            barcodes_path = candidate
            break
    if barcodes_path is None:
        raise FileNotFoundError(f"Missing IsoQuant barcodes file for prefix {prefix}")

    features_path = None
    for suffix in (".features.tsv", ".features.tsv.gz"):
        candidate = Path(f"{prefix}{suffix}")
        if candidate.exists():
            features_path = candidate
            break
    if features_path is None:
        raise FileNotFoundError(f"Missing IsoQuant features file for prefix {prefix}")

    feature_ids = []
    for line in read_lines(features_path):
        parts = line.split("\t")
        feature_ids.append(parts[0])

    barcodes = read_lines(barcodes_path)
    nrows, _, entries = read_mtx(matrix_path)
    col_idx = pick_barcode_index(barcodes, barcode)
    counts = counts_for_column(nrows, entries, col_idx)
    if len(counts) != len(feature_ids):
        raise ValueError(f"Feature/count length mismatch for {prefix}")
    return counts, feature_ids, barcodes, str(prefix)


def find_first(root: Path, patterns: Iterable[str]) -> Optional[Path]:
    if not root.exists():
        return None
    for pattern in patterns:
        hits = sorted(root.glob(pattern))
        if hits:
            return hits[0]
    return None


def resolve_epi2me(outdir: Path, sample: str) -> Tuple[Optional[Path], Optional[Path]]:
    base = outdir / "epi2me" / "create_matrix" / sample
    gene = base / f"{sample}_gene_bc_matrix"
    tx = base / f"{sample}_transcript_bc_matrix"
    if gene.exists() and tx.exists():
        return gene, tx
    gene = find_first(outdir / "epi2me", [f"**/{sample}_gene_bc_matrix", "**/*_gene_bc_matrix"])
    tx = find_first(outdir / "epi2me", [f"**/{sample}_transcript_bc_matrix", "**/*_transcript_bc_matrix"])
    return gene, tx


def resolve_isoquant(outdir: Path, sample: str) -> Tuple[Optional[Path], Optional[Path]]:
    base = outdir / "isoquant" / sample
    gene = find_first(
        base,
        [
            f"**/{sample}.gene_grouped_counts.tsv",
            "**/gene_grouped_counts.tsv",
            "**/gene_grouped_tag_CB_counts.tsv",
            f"**/{sample}.gene_counts.tsv",
        ],
    )
    tx = find_first(
        base,
        [
            f"**/{sample}.transcript_grouped_counts.tsv",
            "**/transcript_grouped_counts.tsv",
            "**/transcript_grouped_tag_CB_counts.tsv",
            f"**/{sample}.transcript_counts.tsv",
        ],
    )
    return gene, tx


def resolve_oarfish(
    outdir: Path, sample: str
) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    base = outdir / "oarfish" / sample
    gene = base / f"{sample}_gene_bc_matrix"
    tx = base / f"{sample}_transcript_bc_matrix"
    quant = base / "quant"
    if not gene.exists():
        gene = find_first(outdir / "oarfish", [f"**/{sample}_gene_bc_matrix", "**/*_gene_bc_matrix"])
    if not tx.exists():
        tx = find_first(
            outdir / "oarfish",
            [f"**/{sample}_transcript_bc_matrix", "**/*_transcript_bc_matrix"],
        )
    if tx is None and quant.exists():
        tx = quant
    return gene, tx, quant


def summarize_method(
    method: str,
    sample: str,
    selected_barcode: str,
    gene_counts: Optional[List[float]],
    tx_counts: Optional[List[float]],
    source_gene: str,
    source_transcript: str,
    n_barcodes: int = 1,
    note: str = "",
) -> MethodMetrics:
    gene_stats = detection_metrics(gene_counts) if gene_counts is not None else None
    tx_stats = detection_metrics(tx_counts) if tx_counts is not None else None
    status = "ok" if gene_stats and tx_stats else "partial" if gene_stats or tx_stats else "missing"
    return MethodMetrics(
        method=method,
        sample=sample,
        barcode=selected_barcode,
        n_barcodes=n_barcodes,
        genes_gt_0=int(gene_stats["detected_gt_0"]) if gene_stats else 0,
        genes_gt_10=int(gene_stats["detected_gt_10"]) if gene_stats else 0,
        transcripts_gt_0=int(tx_stats["detected_gt_0"]) if tx_stats else 0,
        transcripts_gt_10=int(tx_stats["detected_gt_10"]) if tx_stats else 0,
        total_gene_counts=gene_stats["total_counts"] if gene_stats else 0.0,
        total_transcript_counts=tx_stats["total_counts"] if tx_stats else 0.0,
        median_gene_count_detected=gene_stats["median_detected"] if gene_stats else 0.0,
        median_transcript_count_detected=tx_stats["median_detected"] if tx_stats else 0.0,
        source_gene=source_gene,
        source_transcript=source_transcript,
        status=status,
        note=note,
    )


def load_epi2me(
    outdir: Path,
    sample: str,
    barcode: Optional[str],
    gene_dir: Optional[Path] = None,
    tx_dir: Optional[Path] = None,
) -> MethodMetrics:
    if gene_dir is None or tx_dir is None:
        auto_gene, auto_tx = resolve_epi2me(outdir, sample)
        gene_dir = gene_dir or auto_gene
        tx_dir = tx_dir or auto_tx
    if gene_dir is None or tx_dir is None:
        return summarize_method(
            "epi2me", sample, "", None, None, "", "",
            note="Missing create_matrix MEX directories under epi2me/",
        )
    gene_counts, _, gene_src = load_mex_counts(gene_dir, barcode)
    tx_counts, _, tx_src = load_mex_counts(tx_dir, barcode)
    barcodes = read_lines(gene_dir / "barcodes.tsv.gz")
    col_idx = pick_barcode_index(barcodes, barcode)
    return summarize_method(
        "epi2me", sample, barcodes[col_idx], gene_counts, tx_counts,
        gene_src, tx_src, n_barcodes=len(barcodes),
    )


def load_isoquant(outdir: Path, sample: str, barcode: Optional[str]) -> MethodMetrics:
    gene_prefix = resolve_isoquant_mtx_prefix(outdir, sample, "gene")
    tx_prefix = resolve_isoquant_mtx_prefix(outdir, sample, "transcript")
    if gene_prefix is not None and tx_prefix is not None:
        gene_counts, _, barcodes, gene_src = load_isoquant_mtx_counts(gene_prefix, barcode)
        tx_counts, _, _, tx_src = load_isoquant_mtx_counts(tx_prefix, barcode)
        col_idx = pick_barcode_index(barcodes, barcode)
        return summarize_method(
            "isoquant", sample, barcodes[col_idx], gene_counts, tx_counts,
            gene_src, tx_src, n_barcodes=len(barcodes),
        )

    gene_path, tx_path = resolve_isoquant(outdir, sample)
    if gene_path is None or tx_path is None:
        return summarize_method(
            "isoquant", sample, "", None, None, "", "",
            note="Missing grouped MTX or counts TSV under isoquant/",
        )
    gene_counts, gene_src = load_isoquant_grouped(gene_path, barcode)
    tx_counts, tx_src = load_isoquant_grouped(tx_path, barcode)
    with open(gene_path, "rt") as handle:
        header = next(handle).rstrip("\n").split("\t")
    if len(header) == 2 and header[0] == "feature_id":
        selected = barcode or ""
        n_barcodes = 1
        note = "Using per-sample gene_counts.tsv (grouped MTX not found)"
    else:
        barcode_cols = header[1:]
        col_idx = pick_barcode_index(barcode_cols, barcode) if len(barcode_cols) > 1 else 0
        selected = barcode_cols[col_idx]
        n_barcodes = len(barcode_cols)
        note = ""
    return summarize_method(
        "isoquant", sample, selected, gene_counts, tx_counts,
        gene_src, tx_src, n_barcodes=n_barcodes, note=note,
    )


def load_oarfish(
    outdir: Path,
    sample: str,
    barcode: Optional[str],
    gtf: Optional[Path],
    gene_dir: Optional[Path] = None,
    tx_dir: Optional[Path] = None,
) -> MethodMetrics:
    auto_gene, auto_tx, quant_dir = resolve_oarfish(outdir, sample)
    gene_dir = gene_dir or auto_gene
    tx_dir = tx_dir or auto_tx
    note = ""

    if tx_dir is None:
        return summarize_method(
            "oarfish", sample, "", None, None, "", "",
            note="Missing transcript MEX or quant/ directory under oarfish/",
        )

    tx_counts, tx_ids, tx_src = load_mex_counts(tx_dir, barcode)
    barcodes = read_lines(tx_dir / "barcodes.tsv.gz")
    col_idx = pick_barcode_index(barcodes, barcode)
    selected = barcodes[col_idx]

    gene_counts: Optional[List[float]] = None
    gene_src = ""
    if gene_dir is not None and gene_dir.exists():
        gene_counts, _, gene_src = load_mex_counts(gene_dir, barcode)
    elif gtf is not None and gtf.exists():
        tx_to_gene = parse_gtf_transcript_gene(gtf)
        gene_counts = aggregate_transcript_to_gene(tx_counts, tx_ids, tx_to_gene)
        gene_src = f"{tx_src} + GTF aggregation ({gtf})"
        note = "Gene counts aggregated from transcript quant (gene_bc_matrix not published)"
    else:
        note = "Gene metrics unavailable (no gene_bc_matrix; pass --gtf to aggregate from transcripts)"

    return summarize_method(
        "oarfish", sample, selected, gene_counts, tx_counts,
        gene_src, tx_src, n_barcodes=len(barcodes), note=note,
    )


def write_tables(rows: List[MethodMetrics], output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    long_path = output_prefix.with_name(output_prefix.name + "_long.tsv")
    wide_path = output_prefix.with_name(output_prefix.name + "_wide.tsv")

    metric_fields = [
        "genes_gt_0", "genes_gt_10", "transcripts_gt_0", "transcripts_gt_10",
        "total_gene_counts", "total_transcript_counts",
        "median_gene_count_detected", "median_transcript_count_detected",
        "n_barcodes", "barcode", "status", "source_gene", "source_transcript", "note",
    ]

    with open(long_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["method", "sample", "metric", "value"])
        for row in rows:
            data = asdict(row)
            for field in metric_fields:
                writer.writerow([row.method, row.sample, field, data[field]])

    by_method = {row.method: row for row in rows}
    with open(wide_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric"] + list(METHODS))
        for field in metric_fields:
            writer.writerow([field] + [getattr(by_method.get(m), field, "") for m in METHODS])

    print(f"Wrote {long_path}", file=sys.stderr)
    print(f"Wrote {wide_path}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Default pipeline output directory for all methods (if method-specific outdirs omitted)",
    )
    parser.add_argument(
        "--epi2me-outdir",
        type=Path,
        default=None,
        help="Epi2me run outdir (default: --outdir)",
    )
    parser.add_argument(
        "--isoquant-outdir",
        type=Path,
        default=None,
        help="IsoQuant run outdir (default: --outdir)",
    )
    parser.add_argument(
        "--oarfish-outdir",
        type=Path,
        default=None,
        help="Oarfish run outdir (default: --outdir)",
    )
    parser.add_argument("--sample", default="barcode05", help="Sample id (default: barcode05)")
    parser.add_argument(
        "--barcode",
        default=None,
        help="Cell barcode to summarize (optional if matrix has one barcode)",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=None,
        help="Reference GTF for oarfish gene aggregation when gene_bc_matrix is absent",
    )
    parser.add_argument(
        "--epi2me-gene-matrix-dir",
        type=Path,
        default=None,
        help="Override path to epi2me {sample}_gene_bc_matrix directory",
    )
    parser.add_argument(
        "--epi2me-transcript-matrix-dir",
        type=Path,
        default=None,
        help="Override path to epi2me {sample}_transcript_bc_matrix directory",
    )
    parser.add_argument(
        "--oarfish-gene-matrix-dir",
        type=Path,
        default=None,
        help="Override path to oarfish {sample}_gene_bc_matrix directory",
    )
    parser.add_argument(
        "--oarfish-transcript-matrix-dir",
        type=Path,
        default=None,
        help="Override path to oarfish transcript MEX or quant directory",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output prefix without suffix (default: {outdir}/comparison/{sample}_quant_comparison)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.outdir is None and not any(
        (args.epi2me_outdir, args.isoquant_outdir, args.oarfish_outdir)
    ):
        raise SystemExit("Provide --outdir or at least one of --epi2me-outdir, --isoquant-outdir, --oarfish-outdir")

    epi2me_outdir = args.epi2me_outdir or args.outdir
    isoquant_outdir = args.isoquant_outdir or args.outdir
    oarfish_outdir = args.oarfish_outdir or args.outdir
    comparison_base = args.outdir or args.epi2me_outdir or args.isoquant_outdir or args.oarfish_outdir
    output_prefix = args.output or (comparison_base / "comparison" / f"{args.sample}_quant_comparison")

    rows = [
        load_epi2me(
            epi2me_outdir, args.sample, args.barcode,
            gene_dir=args.epi2me_gene_matrix_dir,
            tx_dir=args.epi2me_transcript_matrix_dir,
        ),
        load_isoquant(isoquant_outdir, args.sample, args.barcode),
        load_oarfish(
            oarfish_outdir, args.sample, args.barcode, args.gtf,
            gene_dir=args.oarfish_gene_matrix_dir,
            tx_dir=args.oarfish_transcript_matrix_dir,
        ),
    ]
    write_tables(rows, output_prefix)

    summary_path = output_prefix.with_name(output_prefix.name + "_summary.json")
    with open(summary_path, "wt") as handle:
        json.dump([asdict(r) for r in rows], handle, indent=2)
    print(f"Wrote {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
