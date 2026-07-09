#!/usr/bin/env python3
"""Split oriented ONT reads into R1 (barcode) and R2 (cDNA) using per-read positions.

Reads the bc_tags.tsv.gz produced by extract_barcode.py (columns:
  read_id, cdna_start, barcode_start, CR, CY, anchor_ed, UR, UY)
and for each read present in the tags file:

  R1 = seq[barcode_start : barcode_start + barcode_length]   (36 bp spatial barcode)
  R2 = seq[cdna_start    : barcode_start]                    (clean cDNA, no adapters)

Reads are dropped when:
  - read_id not in bc_tags (barcode extraction failed)
  - R2 length < min_cdna_len (cDNA too short for reliable alignment)
  - R1 length < barcode_length (truncated read)

The synchronized R1/R2 pairs are written to {prefix}_R1.fastq.gz and
{prefix}_R2.fastq.gz so that downstream QUIK processes them as paired-end input
with barcode_start=0.

UMI tags are embedded in the R2 FASTQ header comment as SAM-style optional fields:
  @read_id UR:Z:<raw_umi> UY:Z:<umi_qual>
When minimap2 is run with -y, these are carried through into the BAM as optional
fields, making them directly available to TAG_BAM without a separate TSV join.
If UMI extraction failed (UR is empty), the tags are omitted from the comment.

Usage:
  split_reads.py oriented.fastq.gz bc_tags.tsv.gz \\
      --barcode_length 36 \\
      --min_cdna_len 50 \\
      --prefix sample1
"""
import argparse
import gzip
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [split_reads] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("split_reads")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", type=Path,
                   help="Oriented input FASTQ (gzipped or plain).")
    p.add_argument("bc_tags", type=Path,
                   help="bc_tags.tsv.gz produced by extract_barcode.py.")
    p.add_argument("--barcode_length", type=int, default=36,
                   help="Spatial Barcode length in bp (default: 36).")
    p.add_argument("--min_cdna_len", type=int, default=50,
                   help="Minimum R2 (cDNA) length in bp to keep a read (default: 50).")
    p.add_argument("--prefix", type=str, default="sample",
                   help="Output filename prefix (default: sample).")
    p.add_argument("--r1_out", type=Path, default=None,
                   help="R1 output path. Defaults to {prefix}_R1.fastq.gz.")
    p.add_argument("--r2_out", type=Path, default=None,
                   help="R2 output path. Defaults to {prefix}_R2.fastq.gz.")
    return p.parse_args()


def load_bc_tags(path: Path) -> dict:
    """Load bc_tags.tsv.gz; return {read_id: (cdna_start, barcode_start, UR, UY)}."""
    tags = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx_id = header.index("read_id")
            idx_cdna = header.index("cdna_start")
            idx_bc = header.index("barcode_start")
        except ValueError as e:
            raise ValueError(f"bc_tags missing expected column: {e}") from e
        idx_ur = header.index("UR") if "UR" in header else None
        idx_uy = header.index("UY") if "UY" in header else None
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_id, idx_cdna, idx_bc):
                continue
            read_id = parts[idx_id]
            try:
                cdna_start = int(parts[idx_cdna])
                barcode_start = int(parts[idx_bc])
            except ValueError:
                continue
            ur = parts[idx_ur] if idx_ur is not None and idx_ur < len(parts) else ""
            uy = parts[idx_uy] if idx_uy is not None and idx_uy < len(parts) else ""
            tags[read_id] = (cdna_start, barcode_start, ur, uy)
    logger.info("Loaded %d read tags from %s", len(tags), path)
    return tags


def iter_fastq(path: Path):
    """Yield (name, header, seq, qual) from gzipped or plain FASTQ."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        while True:
            header = fh.readline().rstrip("\n")
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # '+'
            qual = fh.readline().rstrip("\n")
            if not seq or not qual:
                break
            if not header.startswith("@"):
                continue
            name = header[1:].split()[0]
            yield name, header, seq, qual


def main():
    args = parse_args()
    r1_path = args.r1_out or Path(f"{args.prefix}_R1.fastq.gz")
    r2_path = args.r2_out or Path(f"{args.prefix}_R2.fastq.gz")

    tags = load_bc_tags(args.bc_tags)

    n_total = n_written = n_no_tag = n_short_r2 = n_short_r1 = 0

    with gzip.open(r1_path, "wt") as r1_fh, gzip.open(r2_path, "wt") as r2_fh:
        for name, header, seq, qual in iter_fastq(args.fastq):
            n_total += 1

            if name not in tags:
                n_no_tag += 1
                continue

            cdna_start, barcode_start, ur, uy = tags[name]

            r1_seq = seq[barcode_start: barcode_start + args.barcode_length]
            r1_qual = qual[barcode_start: barcode_start + args.barcode_length]
            if len(r1_seq) < args.barcode_length:
                n_short_r1 += 1
                continue

            # cDNA: from end-of-TSO to start-of-barcode
            # If cdna_start == 0 (UMI extraction failed), fall back to 0
            r2_start = cdna_start if cdna_start > 0 else 0
            r2_seq = seq[r2_start: barcode_start]
            r2_qual = qual[r2_start: barcode_start]
            if len(r2_seq) < args.min_cdna_len:
                n_short_r2 += 1
                continue

            # Strip ONT comment fields from the header, keeping only the bare read
            # UUID as the FASTQ read name.  ONT headers use tab-separated fields
            # after the UUID; QUIK appends its _calledidx_N_<barcode> suffix to
            # whatever it finds as the read name.  Without stripping, the suffix
            # lands inside an internal ONT field (e.g. DS:Z:) instead of on the
            # UUID, making BAM QNAME parsing unreliable.
            #
            # R1: bare "@UUID"  → QUIK produces "@UUID_calledidx_N_<CB>"
            # R2: "@UUID UR:Z:<umi>\tUY:Z:<qual>"
            #     minimap2 -y copies the comment tags into the BAM as optional
            #     fields, while QNAME = UUID_calledidx_N_<CB> (from QUIK R2 mod).
            bare_name = "@" + name   # name is already the UUID (before first whitespace)
            r1_header = bare_name
            if ur:
                r2_header = f"{bare_name} UR:Z:{ur}\tUY:Z:{uy}"
            else:
                r2_header = bare_name

            r1_fh.write(f"{r1_header}\n{r1_seq}\n+\n{r1_qual}\n")
            r2_fh.write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")
            n_written += 1

    logger.info(
        "Done: total=%d  written=%d  no_tag=%d  short_r1=%d  short_r2=%d",
        n_total, n_written, n_no_tag, n_short_r1, n_short_r2,
    )
    if n_written == 0:
        logger.warning("No reads written — check bc_tags and window parameters.")
        sys.exit(1)


if __name__ == "__main__":
    main()
