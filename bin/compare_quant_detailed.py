#!/usr/bin/env python3
"""Extended quant method comparison: mapping funnels and gene-set overlap.

Builds on compare_quant_methods.py. Requires --gtf to normalize epi2me gene symbols
to Ensembl gene_id for fair overlap statistics.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Reuse loaders from the base compare script (same directory).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_quant_methods import (  # noqa: E402
    METHODS,
    find_first,
    load_epi2me,
    load_isoquant,
    load_isoquant_grouped,
    load_isoquant_mtx_counts,
    load_mex_counts,
    load_oarfish,
    parse_gtf_transcript_gene,
    resolve_epi2me,
    resolve_isoquant,
    resolve_isoquant_mtx_prefix,
    resolve_oarfish,
)


def parse_gtf_gene_maps(gtf_path: Path) -> Tuple[Dict[str, str], Dict[str, str], Set[str]]:
    """Return gene_name→gene_id, gene_id→gene_name, and all reference gene_ids."""
    name_to_id: Dict[str, str] = {}
    id_to_name: Dict[str, str] = {}
    all_genes: Set[str] = set()
    with open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            if not gid:
                continue
            gene_id = gid.group(1)
            all_genes.add(gene_id)
            if gname:
                name_to_id[gname.group(1)] = gene_id
                id_to_name[gene_id] = gname.group(1)
    return name_to_id, id_to_name, all_genes


def normalize_gene_label(label: str, name_to_id: Dict[str, str]) -> Optional[str]:
    if label.startswith("ENSG"):
        return label.split(".")[0]
    if label in name_to_id:
        return name_to_id[label]
    return None


def counts_dict_from_pairs(labels: List[str], counts: List[float], name_to_id: Dict[str, str]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    unmapped = 0
    for label, val in zip(labels, counts):
        gid = normalize_gene_label(label, name_to_id)
        if gid is None:
            unmapped += 1
            continue
        out[gid] = out.get(gid, 0.0) + float(val)
    return out


def load_gene_counts_dict(
    method: str,
    outdir: Path,
    sample: str,
    barcode: Optional[str],
    name_to_id: Dict[str, str],
) -> Tuple[Dict[str, float], str]:
    if method == "epi2me":
        gene_dir, _ = resolve_epi2me(outdir, sample)
        if gene_dir is None:
            raise FileNotFoundError(f"epi2me gene matrix missing under {outdir}")
        counts, labels, src = load_mex_counts(gene_dir, barcode)
        return counts_dict_from_pairs(labels, counts, name_to_id), src

    if method == "isoquant":
        prefix = resolve_isoquant_mtx_prefix(outdir, sample, "gene")
        if prefix is not None:
            counts, labels, _, src = load_isoquant_mtx_counts(prefix, barcode)
            return counts_dict_from_pairs(labels, counts, name_to_id), src
        gene_path, _ = resolve_isoquant(outdir, sample)
        if gene_path is None:
            raise FileNotFoundError(f"isoquant gene counts missing under {outdir}")
        counts, src = load_isoquant_grouped(gene_path, barcode)
        with open(gene_path, "rt") as handle:
            header = next(handle).rstrip("\n").split("\t")
        if len(header) == 2 and header[0] == "feature_id":
            labels = []
            vals = []
            with open(gene_path, "rt") as handle:
                reader = __import__("csv").reader(handle, delimiter="\t")
                next(reader)
                for row in reader:
                    if len(row) >= 2:
                        labels.append(row[0])
                        vals.append(float(row[1]))
            return counts_dict_from_pairs(labels, vals, name_to_id), src
        raise ValueError(f"Unexpected isoquant grouped format: {gene_path}")

    if method == "oarfish":
        gene_dir, _, _ = resolve_oarfish(outdir, sample)
        if gene_dir is None or not gene_dir.exists():
            raise FileNotFoundError(f"oarfish gene matrix missing under {outdir}")
        counts, labels, src = load_mex_counts(gene_dir, barcode)
        return counts_dict_from_pairs(labels, counts, name_to_id), src

    raise ValueError(method)


def detected_genes(counts: Dict[str, float], min_count: float = 0.0) -> Set[str]:
    return {gid for gid, val in counts.items() if val > min_count}


def pairwise_set_stats(sets: Dict[str, Set[str]]) -> Dict[str, dict]:
    methods = list(sets)
    result: Dict[str, dict] = {}
    union_all = set().union(*sets.values()) if sets else set()
    result["union_all"] = len(union_all)
    result["intersection_all"] = len(set.intersection(*sets.values())) if sets else 0

    for i, a in enumerate(methods):
        for b in methods[i + 1 :]:
            sa, sb = sets[a], sets[b]
            inter = sa & sb
            uni = sa | sb
            key = f"{a}_vs_{b}"
            result[key] = {
                "intersection": len(inter),
                "union": len(uni),
                "jaccard": len(inter) / len(uni) if uni else 0.0,
                f"{a}_only": len(sa - sb),
                f"{b}_only": len(sb - sa),
                f"pct_of_{a}": 100.0 * len(inter) / len(sa) if sa else 0.0,
                f"pct_of_{b}": 100.0 * len(inter) / len(sb) if sb else 0.0,
            }
    for m in methods:
        others = set().union(*(sets[o] for o in methods if o != m))
        result[f"{m}_unique_vs_others"] = len(sets[m] - others)
    return result


def rank_correlation(x: List[float], y: List[float]) -> float:
    if len(x) < 2:
        return float("nan")

    def ranks(vals: List[float]) -> List[float]:
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        r = [0.0] * len(vals)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg_rank
            i = j + 1
        return r

    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den_x = math.sqrt(sum((a - mx) ** 2 for a in rx))
    den_y = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (den_x * den_y) if den_x and den_y else float("nan")


def count_correlation(
    counts_a: Dict[str, float],
    counts_b: Dict[str, float],
    min_count: float = 0.0,
) -> dict:
    shared = detected_genes(counts_a, min_count) & detected_genes(counts_b, min_count)
    if not shared:
        return {"n_shared_detected": 0, "spearman": float("nan"), "pearson": float("nan")}
    xs = [counts_a[g] for g in shared]
    ys = [counts_b[g] for g in shared]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    pearson_num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    pearson_den = math.sqrt(sum((x - mx) ** 2 for x in xs)) * math.sqrt(sum((y - my) ** 2 for y in ys))
    return {
        "n_shared_detected": len(shared),
        "spearman": rank_correlation(xs, ys),
        "pearson": pearson_num / pearson_den if pearson_den else float("nan"),
    }


def parse_flagstat(path: Path) -> dict:
    stats = {"path": str(path)}
    if not path.exists():
        stats["status"] = "missing"
        return stats
    if path.suffix == ".bam":
        try:
            proc = subprocess.run(
                ["samtools", "flagstat", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
            text = proc.stdout
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            stats["status"] = "error"
            stats["error"] = str(exc)
            return stats
    else:
        text = path.read_text()
    for line in text.splitlines():
        parts = line.split()
        if not parts or not parts[0].isdigit():
            continue
        if " in total" in line:
            stats["total"] = int(parts[0])
        elif "primary mapped" in line:
            stats["primary_mapped"] = int(parts[0])
            m = re.search(r"\(([0-9.]+)%", line)
            if m:
                stats["primary_mapped_pct"] = float(m.group(1))
        elif line.rstrip().endswith("mapped") and "primary" not in line and "mate mapped" not in line:
            stats["mapped"] = int(parts[0])
            m = re.search(r"\(([0-9.]+)%", line)
            if m:
                stats["mapped_pct"] = float(m.group(1))
        elif parts[1] == "primary" and "mapped" not in line and "duplicates" not in line:
            stats["primary"] = int(parts[0])
        elif parts[1] == "secondary":
            stats["secondary"] = int(parts[0])
        elif parts[1] == "supplementary":
            stats["supplementary"] = int(parts[0])
    stats["status"] = "ok"
    return stats


def parse_dedup_log(path: Path) -> dict:
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    out = {"path": str(path), "status": "ok"}
    for line in path.read_text().splitlines():
        if "Number of reads out:" in line:
            out["reads_out"] = int(line.rsplit(":", 1)[-1].strip())
        if "Total number of positions deduplicated:" in line:
            out["positions_deduplicated"] = int(line.rsplit(":", 1)[-1].strip())
    return out


def parse_isoquant_log(path: Path) -> dict:
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    out = {"path": str(path), "status": "ok"}
    for line in path.read_text().splitlines():
        if " primary:" in line:
            out["primary_alignments"] = int(line.rsplit(":", 1)[-1].strip())
        if "Total assignments used for analysis:" in line:
            val = line.rsplit(":", 1)[-1].strip().split(",")[0].strip()
            out["assigned_reads"] = int(val)
        if "  unique:" in line:
            out["unique_assignments"] = int(line.rsplit(":", 1)[-1].strip())
    return out


def parse_restrander_stats(path: Path) -> dict:
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    data = json.loads(path.read_text())
    total = data.get("stats", {}).get("totalReads")
    return {"path": str(path), "status": "ok", "total_reads": total}


def count_fastq_reads(path: Path) -> dict:
    """Count reads in a (possibly gzipped) FASTQ file."""
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    try:
        if path.suffix == ".gz":
            proc = subprocess.run(
                ["bash", "-c", f"gzip -dc {path} | wc -l"],
                check=True,
                capture_output=True,
                text=True,
            )
        else:
            proc = subprocess.run(
                ["wc", "-l", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
        n_lines = int(proc.stdout.strip().split()[0])
        return {"path": str(path), "status": "ok", "total_reads": n_lines // 4}
    except (FileNotFoundError, subprocess.CalledProcessError, OSError, ValueError) as exc:
        return {"status": "error", "path": str(path), "error": str(exc)}


def parse_quik_starsolo(outdir: Path, sample: str) -> dict:
    """Return filtered R2 read count from QUIK_STARSOLO (downstream input)."""
    r2 = outdir / "quik_starsolo" / sample / f"{sample}_R2_filtered.fastq.gz"
    stats = count_fastq_reads(r2)
    stats["r2_path"] = str(r2)
    return stats


def parse_matrix_stats(path: Path) -> dict:
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    data = json.loads(path.read_text())
    data["path"] = str(path)
    data["status"] = "ok"
    return data


def collect_mapping_funnel(
    epi2me_outdir: Path,
    isoquant_outdir: Path,
    oarfish_outdir: Path,
    sample: str,
) -> dict:
    restrander = parse_restrander_stats(
        epi2me_outdir / "restrander" / sample / f"{sample}.stats.json"
    )
    quik = parse_quik_starsolo(epi2me_outdir, sample)
    funnel = {
        "restrander": restrander,
        "quik_starsolo": quik,
        "shared_genome_path": {},
        "epi2me": {},
        "isoquant": {},
        "oarfish": {},
    }

    genome_bam = epi2me_outdir / "minimap2" / sample / f"{sample}.bam"
    funnel["shared_genome_path"]["genome_align"] = parse_flagstat(genome_bam)
    funnel["shared_genome_path"]["tagged_genome"] = parse_flagstat(
        epi2me_outdir / "tag_bam" / sample / f"{sample}.tagged.bam"
    )

    funnel["epi2me"]["primary_filter"] = parse_flagstat(
        epi2me_outdir / "epi2me" / "filter_primary_bam" / sample / f"{sample}.primary.bam"
    )
    funnel["epi2me"]["txome_align"] = parse_flagstat(
        epi2me_outdir / "epi2me" / "minimap2_transcriptome" / sample / f"{sample}.bam"
    )
    funnel["epi2me"]["matrix_stats"] = parse_matrix_stats(
        epi2me_outdir / "epi2me" / "create_matrix" / sample / f"{sample}.matrix_stats.json"
    )

    funnel["isoquant"]["primary_filter"] = parse_flagstat(
        isoquant_outdir / "isoquant" / "filter_primary_bam" / sample / f"{sample}.primary.bam"
    )
    funnel["isoquant"]["dedup"] = parse_dedup_log(
        isoquant_outdir / "isoquant" / "umitools_dedup" / sample / f"{sample}.dedup.log"
    )
    funnel["isoquant"]["dedup_bam"] = parse_flagstat(
        isoquant_outdir / "isoquant" / "umitools_dedup" / sample / f"{sample}.dedup.bam"
    )
    funnel["isoquant"]["isoquant_log"] = parse_isoquant_log(
        isoquant_outdir / "isoquant" / sample / sample / "isoquant.log"
    )

    funnel["oarfish"]["txome_align"] = parse_flagstat(
        oarfish_outdir / "oarfish" / sample / "transcriptome_align" / f"{sample}.bam"
    )
    funnel["oarfish"]["txome_mapped"] = parse_flagstat(
        oarfish_outdir / "oarfish" / sample / "transcriptome_align" / f"{sample}.mapped.bam"
    )
    funnel["oarfish"]["dedup"] = parse_dedup_log(
        oarfish_outdir / "oarfish" / "umitools_dedup" / sample / f"{sample}.dedup.log"
    )
    funnel["oarfish"]["dedup_bam"] = parse_flagstat(
        oarfish_outdir / "oarfish" / "umitools_dedup" / sample / f"{sample}.dedup.bam"
    )

    return funnel


def pct(n: Optional[int], d: Optional[int]) -> str:
    if n is None or d is None or d == 0:
        return "—"
    return f"{100.0 * n / d:.1f}%"


METHOD_COLORS = {
    "epi2me": "#059669",
    "isoquant": "#d97706",
    "oarfish": "#7c3aed",
}


def write_plots(
    report: dict,
    gene_counts: Dict[str, Dict[str, float]],
    corr_gt0: Dict[str, dict],
    output_prefix: Path,
    plots_dir: Path,
) -> List[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plots_dir.mkdir(parents=True, exist_ok=True)
    sample = report["sample"]
    written: List[str] = []

    pairs = [("epi2me", "isoquant"), ("epi2me", "oarfish"), ("isoquant", "oarfish")]
    for a, b in pairs:
        shared = detected_genes(gene_counts[a], 0) & detected_genes(gene_counts[b], 0)
        if not shared:
            continue
        xs = np.array([gene_counts[a][g] for g in shared])
        ys = np.array([gene_counts[b][g] for g in shared])
        corr = corr_gt0[f"{a}_vs_{b}"]

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(np.log1p(xs), np.log1p(ys), s=8, alpha=0.35, c="#334155", edgecolors="none")
        lim = max(np.log1p(xs).max(), np.log1p(ys).max()) * 1.05
        ax.plot([0, lim], [0, lim], "--", color="#94a3b8", linewidth=1, label="y = x")
        ax.set_xlabel(f"{a} gene count (log1p)")
        ax.set_ylabel(f"{b} gene count (log1p)")
        ax.set_title(f"{sample}: shared genes (n={len(shared):,})")
        ax.text(
            0.03, 0.97,
            f"Spearman ρ = {corr['spearman']:.3f}\nPearson r = {corr['pearson']:.3f}",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#cbd5e1"),
        )
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal")
        fig.tight_layout()
        out = plots_dir / f"{sample}_scatter_{a}_vs_{b}.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        written.append(str(out))

    by_method = {r["method"]: r for r in report["basic_metrics"]}
    metrics = [
        ("genes_gt_0", "Genes (>0)"),
        ("genes_gt_10", "Genes (>10)"),
        ("transcripts_gt_0", "Transcripts (>0)"),
    ]
    x = np.arange(len(metrics))
    width = 0.25
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, method in enumerate(METHODS):
        vals = [by_method[method][key] for key, _ in metrics]
        ax.bar(x + (i - 1) * width, vals, width, label=method, color=METHOD_COLORS[method])
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in metrics])
    ax.set_ylabel("Features detected")
    ax.set_title(f"{sample}: detection by method")
    ax.legend()
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    out = plots_dir / f"{sample}_bar_detection.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    written.append(str(out))

    overlap = report["gene_overlap"]["gt_0"]
    overlap_labels = [
        "All three",
        "epi2me only",
        "isoquant only",
        "oarfish only",
    ]
    overlap_vals = [
        overlap["intersection_all"],
        overlap.get("epi2me_unique_vs_others", 0),
        overlap.get("isoquant_unique_vs_others", 0),
        overlap.get("oarfish_unique_vs_others", 0),
    ]
    overlap_colors = ["#2563eb", METHOD_COLORS["epi2me"], METHOD_COLORS["isoquant"], METHOD_COLORS["oarfish"]]
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(overlap_labels, overlap_vals, color=overlap_colors)
    ax.set_ylabel("Genes")
    ax.set_title(f"{sample}: gene set overlap (count > 0)")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    for bar, val in zip(bars, overlap_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{val:,}",
                ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    out = plots_dir / f"{sample}_bar_gene_overlap.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    written.append(str(out))

    funnel = report["mapping_funnel"]
    restrander = funnel.get("restrander", {}).get("total_reads")
    quik_reads = funnel.get("quik_starsolo", {}).get("total_reads")
    funnel_rows = [
        ("Restrander", restrander, "#64748b", "restrander"),
        ("QUIK_STARSOLO (R2 filtered)", quik_reads, "#475569", "quik"),
        ("Genome primary mapped", funnel["shared_genome_path"].get("genome_align", {}).get("primary_mapped"), "#2563eb", "downstream"),
        ("Primary filter", funnel["epi2me"].get("primary_filter", {}).get("total"), "#2563eb", "downstream"),
        ("IsoQuant dedup", funnel["isoquant"].get("dedup", {}).get("reads_out"), METHOD_COLORS["isoquant"], "downstream"),
        ("IsoQuant assigned", funnel["isoquant"].get("isoquant_log", {}).get("assigned_reads"), METHOD_COLORS["isoquant"], "downstream"),
        ("Oarfish txome primary", funnel["oarfish"].get("txome_mapped", {}).get("primary_mapped"), METHOD_COLORS["oarfish"], "downstream"),
        ("Oarfish dedup primary", funnel["oarfish"].get("dedup_bam", {}).get("primary_mapped"), METHOD_COLORS["oarfish"], "downstream"),
        ("epi2me gene-tagged", funnel["epi2me"].get("matrix_stats", {}).get("gene_tagged"), METHOD_COLORS["epi2me"], "downstream"),
    ]
    rows = [(l, v, c, tier) for l, v, c, tier in funnel_rows if v is not None]
    labels, vals, colors, tiers = zip(*rows)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ypos = np.arange(len(labels))
    ax.barh(ypos, vals, color=colors)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Reads")
    ax.set_title(f"{sample}: mapping / quantification funnel")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    for i, (v, tier) in enumerate(zip(vals, tiers)):
        if tier == "quik" and restrander:
            ax.text(v, i, f"  {100.0 * v / restrander:.1f}% of restrander", va="center", fontsize=8, color="#475569")
        elif tier == "downstream" and quik_reads:
            ax.text(v, i, f"  {100.0 * v / quik_reads:.1f}% of QUIK", va="center", fontsize=8, color="#475569")
    fig.tight_layout()
    out = plots_dir / f"{sample}_bar_mapping_funnel.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    written.append(str(out))

    return written


def write_markdown(report: dict, path: Path) -> None:
    basic = report["basic_metrics"]
    overlap_gt0 = report["gene_overlap"]["gt_0"]
    overlap_gt10 = report["gene_overlap"]["gt_10"]
    corr_gt0 = report["count_correlation"]["gt_0"]
    funnel = report["mapping_funnel"]
    restrander_reads = funnel.get("restrander", {}).get("total_reads")
    quik_reads = funnel.get("quik_starsolo", {}).get("total_reads")
    out_parent = Path(report["output_prefix"]).parent

    lines = [
        f"# Detailed quantification comparison — `{report['sample']}`",
        "",
        f"*Generated by `compare_quant_detailed.py` · {out_parent}*",
        "",
        "## Detection summary (from basic comparison)",
        "",
        "| Metric | epi2me | isoquant | oarfish |",
        "|--------|--------|----------|---------|",
    ]
    by_method = {r["method"]: r for r in basic}
    for label, key in [
        ("Genes detected (>0)", "genes_gt_0"),
        ("Genes detected (>10)", "genes_gt_10"),
        ("Transcripts detected (>0)", "transcripts_gt_0"),
        ("Total gene counts", "total_gene_counts"),
    ]:
        lines.append(
            "| "
            + label
            + " | "
            + " | ".join(str(by_method[m].get(key, "—")) for m in METHODS)
            + " |"
        )

    lines += [
        "",
        "## Mapping and read funnel",
        "",
        f"Restrander output: **{restrander_reads:,}** reads (shared upstream input).",
        f"QUIK_STARSOLO filtered R2: **{quik_reads:,}** reads ({pct(quik_reads, restrander_reads)} of restrander; input to downstream analysis).",
        "",
        "### Shared genome alignment (epi2me + isoquant)",
        "",
        "| Stage | Reads | % of QUIK | % of restrander | Notes |",
        "|-------|------:|----------:|----------------:|-------|",
    ]
    ga = funnel["shared_genome_path"].get("genome_align", {})
    pf_epi = funnel["epi2me"].get("primary_filter", {})
    lines.append(
        f"| Genome align (all records) | {ga.get('total', '—'):,} | — | — | incl. secondary/supplementary |"
        if ga.get("total") else "| Genome align | — | — | — | |"
    )
    if ga.get("primary_mapped") is not None:
        lines.append(
            f"| Genome primary mapped | {ga['primary_mapped']:,} | {pct(ga['primary_mapped'], quik_reads)} | {pct(ga['primary_mapped'], restrander_reads)} | minimap2 splice · R2 |"
        )
    if pf_epi.get("total") is not None:
        lines.append(
            f"| After primary filter (epi2me/isoquant) | {pf_epi['total']:,} | {pct(pf_epi['total'], quik_reads)} | {pct(pf_epi['total'], restrander_reads)} | -F 0x904 |"
        )

    dedup_iq = funnel["isoquant"].get("dedup", {})
    dedup_iq_bam = funnel["isoquant"].get("dedup_bam", {})
    iq_log = funnel["isoquant"].get("isoquant_log", {})
    lines += [
        "",
        "### isoquant branch",
        "",
        "| Stage | Reads | % of QUIK | % of restrander |",
        "|-------|------:|----------:|----------------:|",
    ]
    if dedup_iq.get("reads_out") is not None:
        lines.append(
            f"| UMI-tools dedup output | {dedup_iq['reads_out']:,} | {pct(dedup_iq['reads_out'], quik_reads)} | {pct(dedup_iq['reads_out'], restrander_reads)} |"
        )
    if iq_log.get("assigned_reads") is not None:
        lines.append(
            f"| IsoQuant assigned reads | {iq_log['assigned_reads']:,} | {pct(iq_log['assigned_reads'], quik_reads)} | {pct(iq_log['assigned_reads'], restrander_reads)} |"
        )

    oar_tx = funnel["oarfish"].get("txome_align", {})
    oar_map = funnel["oarfish"].get("txome_mapped", {})
    oar_dedup = funnel["oarfish"].get("dedup", {})
    oar_dedup_bam = funnel["oarfish"].get("dedup_bam", {})
    lines += [
        "",
        "### oarfish branch (transcriptome path, no genome align)",
        "",
        "| Stage | Primary reads | % of QUIK | % of restrander | Notes |",
        "|-------|-------------:|----------:|----------------:|-------|",
    ]
    if oar_tx.get("primary_mapped") is not None:
        lines.append(
            f"| Txome align primary mapped | {oar_tx['primary_mapped']:,} | {pct(oar_tx['primary_mapped'], quik_reads)} | {pct(oar_tx['primary_mapped'], restrander_reads)} | map-ont · -N 100 |"
        )
    if oar_map.get("primary_mapped") is not None:
        lines.append(
            f"| After mapped filter (-F 4) | {oar_map['primary_mapped']:,} | {pct(oar_map['primary_mapped'], quik_reads)} | {pct(oar_map['primary_mapped'], restrander_reads)} | |"
        )
    if oar_dedup_bam.get("primary_mapped") is not None:
        lines.append(
            f"| UMI-tools dedup (primary) | {oar_dedup_bam['primary_mapped']:,} | {pct(oar_dedup_bam['primary_mapped'], quik_reads)} | {pct(oar_dedup_bam['primary_mapped'], restrander_reads)} | retains secondary records |"
        )

    ms = funnel["epi2me"].get("matrix_stats", {})
    tx_epi = funnel["epi2me"].get("txome_align", {})
    lines += [
        "",
        "### epi2me branch (additional stages)",
        "",
        "| Stage | Reads | Notes |",
        "|-------|------:|-------|",
    ]
    if tx_epi.get("primary_mapped") is not None:
        lines.append(f"| Txome align primary mapped | {tx_epi['primary_mapped']:,} | sample StringTie txome |")
    if ms.get("valid_barcodes") is not None:
        lines.append(f"| Reads with valid barcodes (assign_features) | {ms['valid_barcodes']:,} | from matrix_stats.json |")
    if ms.get("gene_tagged") is not None:
        lines.append(f"| Gene-tagged reads | {ms['gene_tagged']:,} | |")

    lines += [
        "",
        "## Gene detection overlap (normalized to Ensembl `gene_id`)",
        "",
        "epi2me gene symbols mapped via reference GTF `gene_name`. Overlap computed on detected gene sets.",
        "",
        "### Genes with count > 0",
        "",
        f"- **Union (all methods):** {overlap_gt0['union_all']:,}",
        f"- **Intersection (all three):** {overlap_gt0['intersection_all']:,}",
        f"- **Unique to epi2me:** {overlap_gt0.get('epi2me_unique_vs_others', '—'):,}",
        f"- **Unique to isoquant:** {overlap_gt0.get('isoquant_unique_vs_others', '—'):,}",
        f"- **Unique to oarfish:** {overlap_gt0.get('oarfish_unique_vs_others', '—'):,}",
        "",
        "| Pair | Intersection | Union | Jaccard | A only | B only | % of A shared | % of B shared |",
        "|------|-------------:|------:|--------:|-------:|-------:|--------------:|--------------:|",
    ]
    for a, b in [("epi2me", "isoquant"), ("epi2me", "oarfish"), ("isoquant", "oarfish")]:
        key = f"{a}_vs_{b}"
        s = overlap_gt0[key]
        lines.append(
            f"| {a} ∩ {b} | {s['intersection']:,} | {s['union']:,} | {s['jaccard']:.3f} | "
            f"{s[f'{a}_only']:,} | {s[f'{b}_only']:,} | {s[f'pct_of_{a}']:.1f}% | {s[f'pct_of_{b}']:.1f}% |"
        )

    lines += [
        "",
        "### Genes with count > 10",
        "",
        f"- **Intersection (all three):** {overlap_gt10['intersection_all']:,}",
        "",
        "| Pair | Intersection | Jaccard |",
        "|------|-------------:|--------:|",
    ]
    for a, b in [("epi2me", "isoquant"), ("epi2me", "oarfish"), ("isoquant", "oarfish")]:
        key = f"{a}_vs_{b}"
        s = overlap_gt10[key]
        lines.append(f"| {a} ∩ {b} | {s['intersection']:,} | {s['jaccard']:.3f} |")

    lines += [
        "",
        "## Count correlation on shared detected genes (count > 0)",
        "",
        "| Pair | Shared genes | Spearman ρ | Pearson r |",
        "|------|-------------:|-----------:|----------:|",
    ]
    for a, b in [("epi2me", "isoquant"), ("epi2me", "oarfish"), ("isoquant", "oarfish")]:
        c = corr_gt0[f"{a}_vs_{b}"]
        lines.append(
            f"| {a} vs {b} | {c['n_shared_detected']:,} | {c['spearman']:.3f} | {c['pearson']:.3f} |"
        )

    ref = report.get("reference_genes")
    if ref is not None:
        lines += [
            "",
            "## Reference coverage",
            "",
            f"Reference genes in GTF: **{ref['n_reference_genes']:,}**",
            "",
            "| Method | Detected (>0) | % of reference |",
            "|--------|-------------:|---------------:|",
        ]
        for m in METHODS:
            n = overlap_gt0.get(f"{m}_detected", 0)
            lines.append(f"| {m} | {n:,} | {100.0 * n / ref['n_reference_genes']:.1f}% |")

    lines += [
        "",
        "## Plots",
        "",
    ]
    plot_files = report.get("plots", [])
    if plot_files:
        md_dir = path.parent.resolve()
        for plot_path in plot_files:
            plot_name = Path(plot_path).name
            try:
                rel = Path(plot_path).resolve().relative_to(md_dir).as_posix()
            except ValueError:
                rel = plot_name
            lines.append(f"![{plot_name}]({rel})")
            lines.append("")
    else:
        lines.append("_No plots generated._")
        lines.append("")

    lines += [
        "## Interpretation notes",
        "",
        "- **Mapping rates** use QUIK_STARSOLO filtered R2 as the downstream input denominator; restrander is shown for upstream context.",
        "- **Gene overlap** uses Ensembl IDs; epi2me symbols without GTF mapping are excluded from overlap stats.",
        "- **Count correlations** are indicative only — methods use different count types (UMI vs assignment vs EM).",
        "",
        "## Output files",
        "",
        f"- `{path.name}`",
        f"- `{path.with_suffix('.json').name}`",
    ]
    if plot_files:
        lines.append(f"- `{Path(plot_files[0]).parent.name}/` — scatter and bar charts")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epi2me-outdir", type=Path, required=True)
    parser.add_argument("--isoquant-outdir", type=Path, required=True)
    parser.add_argument("--oarfish-outdir", type=Path, required=True)
    parser.add_argument("--sample", default="barcode05")
    parser.add_argument("--barcode", default=None)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=None,
        help="Directory for PNG plots (default: {output_dir}/plots)",
    )
    parser.add_argument("--no-plots", action="store_true", help="Skip plot generation")
    args = parser.parse_args()

    name_to_id, _, ref_genes = parse_gtf_gene_maps(args.gtf)

    basic_rows = [
        load_epi2me(args.epi2me_outdir, args.sample, args.barcode),
        load_isoquant(args.isoquant_outdir, args.sample, args.barcode),
        load_oarfish(args.oarfish_outdir, args.sample, args.barcode, args.gtf),
    ]

    gene_counts: Dict[str, Dict[str, float]] = {}
    sources: Dict[str, str] = {}
    for method, outdir in [
        ("epi2me", args.epi2me_outdir),
        ("isoquant", args.isoquant_outdir),
        ("oarfish", args.oarfish_outdir),
    ]:
        gene_counts[method], sources[method] = load_gene_counts_dict(
            method, outdir, args.sample, args.barcode, name_to_id
        )

    sets_gt0 = {m: detected_genes(gene_counts[m], 0) for m in METHODS}
    sets_gt10 = {m: detected_genes(gene_counts[m], 10) for m in METHODS}

    overlap_gt0 = pairwise_set_stats(sets_gt0)
    overlap_gt10 = pairwise_set_stats(sets_gt10)
    for m in METHODS:
        overlap_gt0[f"{m}_detected"] = len(sets_gt0[m])
        overlap_gt10[f"{m}_detected"] = len(sets_gt10[m])

    corr_gt0 = {}
    for a, b in [("epi2me", "isoquant"), ("epi2me", "oarfish"), ("isoquant", "oarfish")]:
        corr_gt0[f"{a}_vs_{b}"] = count_correlation(gene_counts[a], gene_counts[b], 0)

    funnel = collect_mapping_funnel(
        args.epi2me_outdir, args.isoquant_outdir, args.oarfish_outdir, args.sample
    )

    from dataclasses import asdict

    report = {
        "sample": args.sample,
        "output_prefix": str(args.output),
        "basic_metrics": [asdict(r) for r in basic_rows],
        "gene_count_sources": sources,
        "gene_overlap": {"gt_0": overlap_gt0, "gt_10": overlap_gt10},
        "count_correlation": {"gt_0": corr_gt0},
        "mapping_funnel": funnel,
        "reference_genes": {"n_reference_genes": len(ref_genes)},
        "plots": [],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plots_dir = args.plots_dir or (args.output.parent / "plots")
    if not args.no_plots:
        try:
            report["plots"] = write_plots(report, gene_counts, corr_gt0, args.output, plots_dir)
            for p in report["plots"]:
                print(f"Wrote {p}", file=sys.stderr)
        except ImportError as exc:
            print(f"Skipping plots (matplotlib unavailable): {exc}", file=sys.stderr)

    json_path = args.output.with_suffix(".json")
    md_path = args.output.with_suffix(".md")
    json_path.write_text(json.dumps(report, indent=2))
    write_markdown(report, md_path)
    print(f"Wrote {json_path}", file=sys.stderr)
    print(f"Wrote {md_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
