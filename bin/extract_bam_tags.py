#!/usr/bin/env python3
"""Extract tag TSV files from a tagged genome-alignment BAM.

Outputs (zstd-compressed TSV):
  mapq_tags.tsv.zst     - read_id, mapq  (primary alignments)
  barcode_tags.tsv.zst  - read_id, CR, CY, UR, UY, CB, UB, chr, start, end, SA
                        (CB-sorted for create_matrix chunking)

Usage:
  extract_bam_tags.py tagged.bam --prefix sample1
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [extract_bam_tags] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("extract_bam_tags")

TAGS = ("CR", "CY", "UR", "UY", "CB")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("bam", type=Path, help="Coordinate-sorted tagged BAM.")
    p.add_argument("--prefix", type=str, default="sample",
                   help="Output filename prefix (default: sample).")
    p.add_argument("--threads", type=int, default=2,
                   help="Threads for BAM I/O (default: 2).")
    return p.parse_args()


def get_tag(align, tag):
    """Return BAM tag as string, or '-' if missing."""
    try:
        return align.get_tag(tag)
    except KeyError:
        return "-"


def get_ub(align):
    """Return corrected UMI tag, falling back to raw UMI."""
    ub = get_tag(align, "UB")
    return ub if ub != "-" else get_tag(align, "UR")


def main():
    args = parse_args()
    mapq_rows = []
    barcode_rows = []

    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
        for align in bam.fetch(until_eof=True):
            if align.is_secondary or align.is_supplementary:
                continue
            read_id = align.query_name
            mapq_rows.append((read_id, align.mapping_quality))
            barcode_rows.append({
                "read_id": read_id,
                "CR": get_tag(align, "CR"),
                "CY": get_tag(align, "CY"),
                "UR": get_tag(align, "UR"),
                "UY": get_tag(align, "UY"),
                "CB": get_tag(align, "CB"),
                "UB": get_ub(align),
                "chr": align.reference_name or "-",
                "start": align.reference_start,
                "end": align.reference_end,
                "SA": "-",
            })

    logger.info("Writing %d primary alignments.", len(mapq_rows))
    pd.DataFrame(mapq_rows, columns=["read_id", "mapq"]).to_csv(
        f"{args.prefix}.mapq_tags.tsv.zst",
        sep="\t", index=False, compression="zstd",
    )
    pd.DataFrame(barcode_rows).sort_values(
        ["CB", "read_id"], kind="mergesort",
    ).to_csv(
        f"{args.prefix}.barcode_tags.tsv.zst",
        sep="\t", index=False, compression="zstd",
    )


if __name__ == "__main__":
    main()
