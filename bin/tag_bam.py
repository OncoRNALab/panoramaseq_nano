#!/usr/bin/env python3
"""Add standardised barcode and UMI tags to a minimap2-aligned BAM.

For each alignment record:
  - CB: corrected spatial barcode parsed from QUIK read name suffix
  - CR / CY: raw barcode and quality from bc_tags.tsv.gz (joined by read UUID)
  - UR / UY: raw UMI and quality (bc_tags preferred, else existing BAM tags)
  - UB: corrected UMI placeholder (= UR; UMI-tools dedup updates this later)

QUIK read names follow:
  {uuid}_calledidx_{idx}_{barcode}

Usage:
  tag_bam.py aligned.bam bc_tags.tsv.gz --prefix sample1
"""
import argparse
import gzip
import logging
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [tag_bam] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("tag_bam")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam", type=Path, help="Coordinate-sorted input BAM.")
    p.add_argument("bc_tags", type=Path,
                   help="bc_tags.tsv.gz produced by extract_barcode.py.")
    p.add_argument("--prefix", type=str, default="sample",
                   help="Output filename prefix (default: sample).")
    p.add_argument("--barcode_length", type=int, default=36,
                   help="Expected corrected barcode length from QUIK (default: 36).")
    return p.parse_args()


def load_bc_tags(path: Path) -> Dict[str, Tuple[str, str, str, str]]:
    """Load bc_tags.tsv.gz; return {read_id: (CR, CY, UR, UY)}."""
    tags: Dict[str, Tuple[str, str, str, str]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: header.index(name) for name in ("read_id", "CR", "CY", "UR", "UY")}

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx.values()):
                continue
            tags[parts[idx["read_id"]]] = (
                parts[idx["CR"]],
                parts[idx["CY"]],
                parts[idx["UR"]],
                parts[idx["UY"]],
            )

    logger.info("Loaded %d bc_tags entries from %s", len(tags), path)
    return tags


def parse_quik_qname(qname: str, barcode_length: int) -> Tuple[Optional[str], Optional[str]]:
    """Return (read_uuid, corrected_barcode) from a QUIK-formatted QNAME."""
    marker = "_calledidx_"
    if marker not in qname:
        return qname, None

    read_uuid, rest = qname.split(marker, 1)
    if "_" not in rest:
        return read_uuid, None

    _idx, barcode = rest.split("_", 1)
    if len(barcode) != barcode_length:
        logger.debug("Unexpected barcode length %d for %s", len(barcode), qname)
    return read_uuid, barcode


def get_tag_str(read: pysam.AlignedSegment, tag: str) -> str:
    if not read.has_tag(tag):
        return ""
    value = read.get_tag(tag)
    return value if isinstance(value, str) else str(value)


def tag_bam(in_bam: Path, bc_tags: Dict[str, Tuple[str, str, str, str]],
            out_bam: Path, barcode_length: int) -> Dict[str, int]:
    stats = {
        "total": 0,
        "cb_set": 0,
        "cr_set": 0,
        "ur_set": 0,
        "missing_bc_tags": 0,
        "missing_cb": 0,
    }

    with pysam.AlignmentFile(in_bam, "rb") as infile, \
         pysam.AlignmentFile(out_bam, "wb", template=infile) as outfile:
        for read in infile:
            stats["total"] += 1
            read_uuid, cb = parse_quik_qname(read.query_name, barcode_length)

            cr = cy = ur = uy = ""
            if read_uuid in bc_tags:
                cr, cy, ur, uy = bc_tags[read_uuid]
            else:
                stats["missing_bc_tags"] += 1

            if not ur:
                ur = get_tag_str(read, "UR")
            if not uy:
                uy = get_tag_str(read, "UY")

            if cb:
                read.set_tag("CB", cb, "Z", replace=True)
                stats["cb_set"] += 1
            else:
                stats["missing_cb"] += 1

            if cr:
                read.set_tag("CR", cr, "Z", replace=True)
                stats["cr_set"] += 1
            if cy:
                read.set_tag("CY", cy, "Z", replace=True)
            if ur:
                read.set_tag("UR", ur, "Z", replace=True)
                read.set_tag("UB", ur, "Z", replace=True)
                stats["ur_set"] += 1
            if uy:
                read.set_tag("UY", uy, "Z", replace=True)

            outfile.write(read)

    return stats


def main():
    args = parse_args()
    out_bam = Path(f"{args.prefix}.tagged.bam")

    bc_tags = load_bc_tags(args.bc_tags)
    stats = tag_bam(args.bam, bc_tags, out_bam, args.barcode_length)

    pysam.index(str(out_bam))

    logger.info(
        "Tagged BAM written to %s: total=%d cb=%d cr=%d ur=%d missing_bc_tags=%d missing_cb=%d",
        out_bam,
        stats["total"],
        stats["cb_set"],
        stats["cr_set"],
        stats["ur_set"],
        stats["missing_bc_tags"],
        stats["missing_cb"],
    )


if __name__ == "__main__":
    main()
