#!/usr/bin/env python3
"""Extract Spatial Barcodes and structured UMIs from oriented ONT reads.

Strategy (adapted from wf-single-cell extract_barcode.py):
  After orientation, every read has the structure (5'→3'):
    [adapter_anchor(23bp)] - [VVVV TT VVVV TT VVVV TT VVVV] - [TTT GGG] -
    [cDNA] - [polyA(~18bp)] - [SpatialBarcode(36bp)] - [RC(RT-Adapter)(45bp)] - [OuterBC]

  Two Smith-Waterman alignments are performed per read:

  3' alignment (Spatial Barcode):
    - Query: last <window> bases of oriented read
    - Probe: polyA(18) + N*36 + RC(RT-Adapter)
    - Extracts: barcode_start, CR (raw barcode), CY (barcode quality), anchor_ed

  5' alignment (UMI):
    - Query: first <umi_window> bases of oriented read
    - Probe: adapter_anchor only (e.g. TCTGTTGGTGCTGATATTGCTTT, 23 bp)
    - After the anchor end, UMI is at fixed positions in the next 22 bases:
        [VVVV TT VVVV TT VVVV TT VVVV]  → indices [0-3, 6-9, 12-15, 18-21]
    - Extracts: cdna_start, UR (raw 16-base UMI), UY (UMI quality)

Output TSV columns (bc_tags.tsv.gz):
  read_id, cdna_start, barcode_start, CR, CY, anchor_ed, UR, UY

Usage:
  extract_barcode.py oriented.fastq.gz \\
      --rt_adapter CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG \\
      --adapter_anchor TCTGTTGGTGCTGATATTGC \\
      --barcode_length 36 \\
      --window 150 \\
      --umi_window 80 \\
      --max_anchor_ed 5 \\
      --max_adapter_ed 3 \\
      --min_barcode_qv 15 \\
      --output bc_tags.tsv.gz
"""
import argparse
import gzip
import logging
import statistics
import sys
from pathlib import Path

import editdistance
import parasail

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [extract_barcode] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("extract_barcode")

POLYA_LENGTH = 18


def rc(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", type=Path, help="Oriented input FASTQ (gzipped or plain).")
    # 3' alignment params
    p.add_argument(
        "--rt_adapter",
        default="CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        help="Constant 45 nt cDNA-RT adapter anchor sequence.",
    )
    p.add_argument("--barcode_length", type=int, default=36,
                   help="Spatial Barcode length in bp (default: 36).")
    p.add_argument("--polya_length", type=int, default=POLYA_LENGTH,
                   help="Poly-A tail length in the 3' probe (default: 18).")
    p.add_argument("--window", type=int, default=150,
                   help="Bases from 3' end used as alignment query (default: 150).")
    p.add_argument("--max_anchor_ed", type=int, default=5,
                   help="Max edit distance for RC(RT-adapter) anchor check (default: 5).")
    p.add_argument("--min_barcode_qv", type=int, default=15,
                   help="Minimum Phred quality across all barcode bases (default: 15).")
    p.add_argument("--min_barcode_len", type=int, default=None,
                   help="Minimum accepted barcode length before padding. Barcodes shorter "
                        "than this are dropped; barcodes between min_barcode_len and "
                        "barcode_length are right-padded with G's (quality '!'). "
                        "Default: barcode_length (no padding, strict).")
    # 5' alignment params
    p.add_argument(
        "--adapter_anchor",
        default="TCTGTTGGTGCTGATATTGC",
        help="20 bp constant adapter anchor at the 5' end, preceding the UMI (default: in-house).",
    )
    p.add_argument("--umi_window", type=int, default=80,
                   help="Bases from 5' end used for UMI alignment query (default: 80).")
    p.add_argument("--max_adapter_ed", type=int, default=3,
                   help="Max edit distance for 5' adapter anchor check (default: 3). "
                        "Reads exceeding this still get CR/CY but get empty UR/UY.")
    # Shared SW params
    p.add_argument("--gap_open", type=int, default=2,
                   help="Parasail SW gap-open penalty (default: 2).")
    p.add_argument("--gap_extend", type=int, default=4,
                   help="Parasail SW gap-extend penalty (default: 4).")
    p.add_argument("--match", type=int, default=5,
                   help="Base match score (default: 5).")
    p.add_argument("--mismatch", type=int, default=-1,
                   help="Base mismatch score (default: -1).")
    p.add_argument("--acg_to_n_match", type=int, default=1,
                   help="Score for A/C/G <-> N match (default: 1).")
    p.add_argument("--t_to_n_match", type=int, default=1,
                   help="Score for T <-> N match (default: 1).")
    # Outputs
    p.add_argument("--output", type=Path, default=None,
                   help="Output TSV.gz path. Defaults to stdout.")
    p.add_argument("--consensus_output", type=Path, default=None,
                   help="Optional file to write the median barcode_start integer.")
    return p.parse_args()


def build_matrix(match, mismatch, acg_to_n_match, t_to_n_match):
    """Build ACGTN parasail scoring matrix with custom N-match scores."""
    matrix = parasail.matrix_create("ACGTN", match, mismatch)
    # N <-> A/C/G (and symmetric)
    for i in [4, 9, 14, 20, 21, 22]:
        matrix.pointer[0].matrix[i] = acg_to_n_match
    # N <-> T (and T <-> N)
    matrix.pointer[0].matrix[19] = t_to_n_match
    matrix.pointer[0].matrix[23] = t_to_n_match
    return matrix


class FastqRecord:
    def __init__(self, name: str, sequence: str, quality: str):
        self.name = name
        self.sequence = sequence
        self.quality = quality


def iter_fastq(path: Path):
    """Iterate FASTQ records from gzipped or plain input."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                break
            header = header.rstrip("\n")
            if not header.startswith("@"):
                continue
            name = header[1:].split()[0]
            yield FastqRecord(name=name, sequence=seq.rstrip("\n"), quality=qual.rstrip("\n"))


def extract_from_alignment(p_aln, probe_anchor, barcode_length, query_seq, query_qual,
                           min_barcode_len=None):
    """Parse 3' SW traceback. Returns (barcode, barcode_qual, bc_start_in_query, min_qv, anchor_ed).

    If the extracted barcode is shorter than barcode_length but >= min_barcode_len,
    it is right-padded with G's (quality '!' = Phred 0) to reach barcode_length.
    If min_barcode_len is None it defaults to barcode_length (no padding, strict mode).
    """
    if min_barcode_len is None:
        min_barcode_len = barcode_length

    ref_aln = p_aln.traceback.ref
    qry_aln = p_aln.traceback.query

    bc_start_col = ref_aln.find("N")
    if bc_start_col == -1:
        return None, None, -1, 0, len(probe_anchor)

    # Walk to end of N-block
    n_seen = 0
    bc_end_col = bc_start_col
    for i in range(bc_start_col, len(ref_aln)):
        if ref_aln[i] == "N":
            n_seen += 1
        if n_seen == barcode_length:
            bc_end_col = i + 1
            break

    # Anchor check: query bases after the N-block vs. RC(RT-Adapter)
    anchor_hit = qry_aln[bc_end_col:].replace("-", "")[:len(probe_anchor)]
    anchor_ed = editdistance.eval(anchor_hit, probe_anchor)

    # Barcode start in query window via alignment coordinates
    n_query_in_aln = len(qry_aln) - qry_aln.count("-")
    query_aln_start = p_aln.end_query - n_query_in_aln + 1
    prefix_query_bases = len(qry_aln[:bc_start_col].replace("-", ""))
    bc_start_in_query = query_aln_start + prefix_query_bases

    barcode = query_seq[bc_start_in_query: bc_start_in_query + barcode_length]
    actual_len = len(barcode)

    if actual_len < min_barcode_len:
        return None, None, -1, 0, anchor_ed

    # Pad with G's if shorter than barcode_length but >= min_barcode_len
    if actual_len < barcode_length:
        pad = barcode_length - actual_len
        barcode   = barcode   + "G" * pad
        bc_qual_str = query_qual[bc_start_in_query: bc_start_in_query + actual_len] + "!" * pad
    else:
        bc_qual_str = query_qual[bc_start_in_query: bc_start_in_query + barcode_length]

    min_qv = min(ord(c) - 33 for c in bc_qual_str) if bc_qual_str else 0

    return barcode, bc_qual_str, bc_start_in_query, min_qv, anchor_ed


def extract_umi_from_alignment(p_aln, adapter_anchor, query_seq, query_qual):
    """Parse 5' SW anchor alignment; extract UMI at fixed offset after the anchor.

    The adapter_anchor (e.g. TCTGTTGGTGCTGATATTGCTTT) is the last constant
    sequence before the UMI region.  After the anchor the TSO structure is:

        [VVVV TT VVVV TT VVVV TT VVVV] TTT GGG
         0123 45 6789    ...
         <----  22 bp UMI region  ---->  <6 trailing>

    The 16 UMI bases are the non-T positions:
        indices [0,1,2,3, 6,7,8,9, 12,13,14,15, 18,19,20,21]

    cdna_start = first base of cDNA = anchor_end + 1 + 22 (UMI region) + 6 (TTT+GGG)

    Returns ("", "", 0, large_ed) if the read is too short to extract the UMI.
    """
    # Last aligned query position (0-based, inclusive)
    anchor_end = p_aln.end_query

    # Edit distance of the aligned portion vs the full anchor
    qry_aln = p_aln.traceback.query
    anchor_hit = qry_aln.replace("-", "")[-len(adapter_anchor):]
    adapter_ed = editdistance.eval(anchor_hit, adapter_anchor)

    # Fixed-offset extraction — no N-walking needed
    UMI_REGION_LEN = 22    # VVVV TT VVVV TT VVVV TT VVVV
    UMI_POSITIONS  = [0, 1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21]
    TRAILING_CONST = 6     # TTT + GGG after the UMI region

    umi_start = anchor_end + 1
    umi_region      = query_seq [umi_start: umi_start + UMI_REGION_LEN]
    umi_region_qual = query_qual[umi_start: umi_start + UMI_REGION_LEN]

    if len(umi_region) < UMI_REGION_LEN:
        return "", "", 0, len(adapter_anchor)

    umi_seq  = "".join(umi_region     [i] for i in UMI_POSITIONS)
    umi_qual = "".join(umi_region_qual[i] for i in UMI_POSITIONS)
    cdna_start = umi_start + UMI_REGION_LEN + TRAILING_CONST

    return umi_seq, umi_qual, cdna_start, adapter_ed


def main():
    args = parse_args()

    # --- 3' probe (spatial barcode) ---
    rt_adapter_rc = rc(args.rt_adapter)
    probe_seq = "A" * args.polya_length + "N" * args.barcode_length + rt_adapter_rc
    probe_anchor = rt_adapter_rc
    logger.info("3' probe (%d bp): %s", len(probe_seq), probe_seq)

    # --- 5' probe (UMI): align the anchor only ---
    # After the anchor ends, UMI bases are extracted at fixed positions in the
    # next 22 bases (VVVVTTVVVVTTVVVVTTVVVV).  No N-pattern probe is needed.
    umi_probe_seq = args.adapter_anchor
    logger.info("5' UMI anchor probe (%d bp): %s", len(umi_probe_seq), umi_probe_seq)

    matrix = build_matrix(
        args.match, args.mismatch, args.acg_to_n_match, args.t_to_n_match
    )

    min_barcode_len = args.min_barcode_len if args.min_barcode_len is not None else args.barcode_length
    if min_barcode_len > args.barcode_length:
        raise ValueError(f"--min_barcode_len ({min_barcode_len}) cannot exceed --barcode_length ({args.barcode_length})")
    logger.info("Barcode length: %d  min accepted: %d  padding: %s",
                args.barcode_length, min_barcode_len,
                "G's" if min_barcode_len < args.barcode_length else "disabled")

    barcode_starts = []

    def write_records(out_fh):
        out_fh.write("read_id\tcdna_start\tbarcode_start\tCR\tCY\tanchor_ed\tUR\tUY\n")
        n_total = n_bc_pass = n_umi_pass = 0
        for read in iter_fastq(args.fastq):
            n_total += 1
            seq = read.sequence
            qual = read.quality

            # --- 3' alignment: spatial barcode ---
            if len(seq) < args.window:
                query_seq_3 = seq
                query_qual_3 = qual
                read_offset = 0
            else:
                query_seq_3 = seq[-args.window:]
                query_qual_3 = qual[-args.window:]
                read_offset = len(seq) - args.window

            aln_3 = parasail.sw_trace(
                s1=query_seq_3,
                s2=probe_seq,
                open=args.gap_open,
                extend=args.gap_extend,
                matrix=matrix,
            )

            barcode, bc_qual_str, bc_start_in_query, min_qv, anchor_ed = extract_from_alignment(
                aln_3, probe_anchor, args.barcode_length, query_seq_3, query_qual_3,
                min_barcode_len=min_barcode_len
            )

            if barcode is None or anchor_ed > args.max_anchor_ed or min_qv < args.min_barcode_qv:
                continue

            n_bc_pass += 1
            barcode_start_in_read = read_offset + bc_start_in_query
            barcode_starts.append(barcode_start_in_read)

            # --- 5' alignment: UMI ---
            umi_window_seq = seq[:args.umi_window]
            umi_window_qual = qual[:args.umi_window]

            aln_5 = parasail.sw_trace(
                s1=umi_window_seq,
                s2=umi_probe_seq,
                open=args.gap_open,
                extend=args.gap_extend,
                matrix=matrix,
            )

            umi_seq, umi_qual, cdna_start, adapter_ed = extract_umi_from_alignment(
                aln_5, args.adapter_anchor, umi_window_seq, umi_window_qual
            )

            if adapter_ed > args.max_adapter_ed:
                umi_seq = ""
                umi_qual = ""
                cdna_start = 0
            else:
                n_umi_pass += 1

            out_fh.write(
                f"{read.name}\t{cdna_start}\t{barcode_start_in_read}\t"
                f"{barcode}\t{bc_qual_str}\t{anchor_ed}\t{umi_seq}\t{umi_qual}\n"
            )

        logger.info(
            "Reads: total=%d  barcode_pass=%d  umi_pass=%d",
            n_total, n_bc_pass, n_umi_pass,
        )

    if args.output:
        with gzip.open(args.output, "wt") as gz:
            write_records(gz)
    else:
        write_records(sys.stdout)

    if barcode_starts:
        consensus = int(statistics.median(barcode_starts))
        logger.info("Consensus barcode_start (median): %d", consensus)
        if args.consensus_output:
            args.consensus_output.write_text(str(consensus) + "\n")
        print(consensus)
    else:
        logger.warning("No reads passed the barcode extraction filters.")
        if args.consensus_output:
            args.consensus_output.write_text("0\n")
        print(0)


if __name__ == "__main__":
    main()
