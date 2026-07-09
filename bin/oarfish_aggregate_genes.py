#!/usr/bin/env python3
"""Aggregate Oarfish transcript-level MEX matrices to gene level."""
import argparse
import gzip
import json
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path


def read_lines(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line:
                yield line


def parse_gtf_attribute(field: str, key: str):
    match = re.search(rf'{key} "([^"]+)"', field)
    return match.group(1) if match else None


def load_transcript_gene_map(gtf_path: Path):
    mapping = {}
    with open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            transcript_id = parse_gtf_attribute(parts[8], "transcript_id")
            gene_id = parse_gtf_attribute(parts[8], "gene_id")
            if transcript_id and gene_id:
                mapping[transcript_id] = gene_id
    return mapping


def parse_feature(line: str, tx_to_gene: dict):
    if "\t" in line:
        parts = line.split("\t")
        if len(parts) >= 2:
            return parts[0], parts[1]
    if "|" in line:
        parts = line.split("|")
        if len(parts) >= 2:
            return parts[0], parts[1]
    transcript_id = line.strip()
    if not transcript_id:
        raise ValueError("Empty feature line")
    gene_id = tx_to_gene.get(transcript_id, transcript_id)
    return transcript_id, gene_id


def read_mtx(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                nrows, ncols, nentries = map(int, line.split())
                break
        rows, cols, data = [], [], []
        for _ in range(nentries):
            parts = handle.readline().split()
            rows.append(int(parts[0]) - 1)
            cols.append(int(parts[1]) - 1)
            data.append(float(parts[2]))
    return nrows, ncols, rows, cols, data


def write_mtx(path: Path, nrows: int, ncols: int, rows, cols, data):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        handle.write("%%MatrixMarket matrix coordinate real general\n")
        handle.write("%written_by oarfish_aggregate_genes.py\n")
        handle.write(f"{nrows} {ncols} {len(data)}\n")
        for row, col, val in zip(rows, cols, data):
            handle.write(f"{row + 1} {col + 1} {val}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", default="sample")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tx_dir = args.output_dir / f"{args.prefix}_transcript_bc_matrix"
    gene_dir = args.output_dir / f"{args.prefix}_gene_bc_matrix"
    tx_dir.mkdir(exist_ok=True)
    gene_dir.mkdir(exist_ok=True)

    for src, dst in [
        (args.features, tx_dir / "features.tsv.gz"),
        (args.barcodes, tx_dir / "barcodes.tsv.gz"),
        (args.matrix, tx_dir / "matrix.mtx.gz"),
    ]:
        shutil.copy2(src, dst)

    tx_to_gene = load_transcript_gene_map(args.gtf)
    features = list(read_lines(args.features))
    barcodes = list(read_lines(args.barcodes))
    gene_ids = [parse_feature(line, tx_to_gene)[1] for line in features]

    nrows, ncols, rows, cols, data = read_mtx(args.matrix)
    if nrows != len(features) and ncols == len(features):
        rows, cols = cols, rows
        nrows, ncols = ncols, nrows
    if nrows != len(features):
        sys.exit(
            f"Matrix rows {nrows} do not match {len(features)} transcript features"
        )
    if ncols != len(barcodes):
        sys.exit(
            f"Matrix columns {ncols} do not match {len(barcodes)} barcodes"
        )

    gene_to_idx = {gene: idx for idx, gene in enumerate(sorted(set(gene_ids)))}
    row_gene = [gene_to_idx[gene] for gene in gene_ids]

    gene_data = defaultdict(float)
    for row, col, val in zip(rows, cols, data):
        gene_data[(row_gene[row], col)] += val

    gene_rows, gene_cols, gene_vals = zip(*[
        (g, c, v) for (g, c), v in gene_data.items()
    ]) if gene_data else ([], [], [])

    with gzip.open(gene_dir / "features.tsv.gz", "wt") as handle:
        for gene in sorted(gene_to_idx, key=gene_to_idx.get):
            handle.write(f"{gene}\t{gene}\tGene Expression\n")

    with gzip.open(gene_dir / "barcodes.tsv.gz", "wt") as handle:
        for bc in barcodes:
            handle.write(f"{bc}\n")

    write_mtx(
        gene_dir / "matrix.mtx.gz",
        len(gene_to_idx),
        len(barcodes),
        list(gene_rows),
        list(gene_cols),
        list(gene_vals),
    )

    mapped = sum(
        1 for line in features
        if parse_feature(line, tx_to_gene)[0] in tx_to_gene
    )
    stats = {
        "transcripts": len(features),
        "genes": len(gene_to_idx),
        "barcodes": len(barcodes),
        "transcripts_mapped_to_gtf": mapped,
        "transcripts_unmapped_to_gtf": len(features) - mapped,
    }
    stats_path = args.output_dir / f"{args.prefix}.aggregation_stats.json"
    stats_path.write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
