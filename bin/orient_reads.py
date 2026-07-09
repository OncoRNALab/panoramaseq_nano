#!/usr/bin/env python3
"""DEPRECATED: replaced by Restrander (modules/local/restrander). Not used in the main workflow.

Orient ONT reads so the RT-Adapter (and Spatial Barcode) is at the 3' end.

Strategy (adapted from wf-single-cell adapter_scan_vsearch.py):
  1. Run vsearch --usearch_global against a 4-orientation adapter FASTA.
  2. Parse the hit table: determine which adapter target was matched.
  3. Targets whose name ends with '_r' indicate that the read is in antisense
     orientation (the reverse complement of the adapter was found). These reads
     are reverse-complemented so that the RT-Adapter ends up at the 3' end of
     every output record.
  4. Reads with no vsearch hit are passed through unchanged (unstranded).

Expected adapter FASTA format (4 entries):
  >tso_f          <- TSO forward (5' landmark)
  >tso_r          <- TSO reverse complement
  >rt_adapter_f   <- RT-Adapter forward (3' landmark, sense orientation)
  >rt_adapter_r   <- RT-Adapter reverse complement

Usage:
  # Option 1: provide a ready FASTA
  orient_reads.py reads.fastq.gz \\
      --adapters adapters.fasta \\
      --min_adapter_id 0.7 \\
      --threads 4 \\
      --output oriented.fastq.gz

  # Option 2: provide TSO + RT adapter directly
  orient_reads.py reads.fastq.gz \\
      --tso_sequence "TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTmGmGmG" \\
      --rt_adapter "CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG" \\
      --output oriented.fastq.gz
"""
import argparse
import gzip
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [orient_reads] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("orient_reads")

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rev_cmp(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", type=Path, help="Input FASTQ (gzipped or plain).")
    p.add_argument("--adapters", type=Path,
                   help="FASTA with 4-orientation adapter sequences.")
    p.add_argument("--tso_sequence", type=str, default=None,
                   help="TSO/composite sequence. Supports IUPAC and modified-base "
                        "notation (e.g. V, mG). Used when --adapters is not set.")
    p.add_argument("--rt_adapter", type=str,
                   default="CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
                   help="RT-adapter anchor sequence used when generating adapters.")
    p.add_argument("--min_adapter_id", type=float, default=0.7,
                   help="Minimum alignment identity for vsearch (default: 0.7).")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for vsearch (default: 4).")
    p.add_argument("--output", type=Path, default=None,
                   help="Output FASTQ.gz path. Defaults to stdout.")
    p.add_argument("--summary", type=Path, default=None,
                   help="Optional TSV of per-read orientation decisions.")
    return p.parse_args()


def normalize_probe_sequence(seq: str) -> str:
    """Normalize oligo notation to DNA alphabet for sequence search.

    - remove modified-base marks like 'mG' -> 'G'
    - convert U -> T
    - convert IUPAC degenerate codes to N (e.g. V -> N)
    """
    s = seq.strip().upper()
    # remove common modified-base prefix marker "m"
    s = s.replace("M", "")
    s = s.replace("U", "T")
    allowed = set("ACGT")
    return "".join(ch if ch in allowed else "N" for ch in s)


def write_adapters_fasta(tso: str, rt_adapter: str, out_fasta: Path) -> None:
    """Write 4-orientation adapters FASTA from TSO and RT-adapter sequences."""
    tso_n = normalize_probe_sequence(tso)
    rt_n = normalize_probe_sequence(rt_adapter)
    # Naming convention: _f = found in correctly oriented (top-strand) reads → KEEP
    #                    _r = found in bottom-strand reads → RC
    # TSO is at the 5' end of the top strand:
    #   tso_f  = TSO forward sequence  → top-strand read → KEEP
    #   tso_r  = RC(TSO)               → bottom-strand read → RC
    # RT-adapter is at the 3' end of the top strand as RC(rt_adapter):
    #   rt_adapter_f = RC(rt_adapter)  → found at 3' of top-strand read → KEEP
    #   rt_adapter_r = rt_adapter_f    → found at 5' of bottom-strand read → RC
    adapters = {
        "tso_f": tso_n,
        "tso_r": rev_cmp(tso_n),
        "rt_adapter_f": rev_cmp(rt_n),
        "rt_adapter_r": rt_n,
    }
    with open(out_fasta, "w") as fh:
        for name, seq in adapters.items():
            fh.write(f">{name}\n{seq}\n")
    logger.info("Generated adapters FASTA at %s", out_fasta)


def run_vsearch(fastq: Path, adapters: Path, min_id: float,
                threads: int, out_tsv: Path) -> None:
    """Run vsearch usearch_global and write hits to out_tsv."""
    fields = "query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl"
    cmd = (
        f"seqkit fq2fa {fastq} | vsearch"
        f" --usearch_global -"
        f" --db {adapters}"
        f" --minseqlength 20"
        f" --maxaccepts 5"
        f" --id {min_id}"
        f" --strand plus"
        f" --wordlength 3"
        f" --minwordmatches 10"
        f" --output_no_hits"
        f" --userfields '{fields}'"
        f" --userout {out_tsv}"
        f" --threads {threads}"
    )
    logger.info("Running vsearch: %s", cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        logger.error("vsearch stderr:\n%s", result.stderr.decode())
        raise RuntimeError("vsearch failed")
    logger.info("vsearch finished.")


def parse_hits(tsv: Path) -> dict:
    """Return {read_id: best_target} for all reads."""
    cols = [
        "query", "target", "id", "alnlen", "mism", "opens",
        "qilo", "qihi", "qstrand", "tilo", "tihi", "ql", "tl",
    ]
    best = {}  # read_id -> best (id_score, target)
    with open(tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(cols, parts))
            read_id = row["query"]
            if row["target"] == "*":
                # no hit
                if read_id not in best:
                    best[read_id] = (0.0, "*")
                continue
            score = float(row["id"])
            if read_id not in best or score > best[read_id][0]:
                best[read_id] = (score, row["target"])
    return {k: v[1] for k, v in best.items()}


def needs_rc(target: str) -> bool:
    """Return True if the read should be reverse-complemented.

    Targets ending with '_r' mean the reverse-complement adapter was found on the
    read → the read is in antisense orientation → RC to put RT-Adapter at 3' end.
    """
    return target.endswith("_r")


def open_fastq(path: Path):
    """Open a (possibly gzipped) FASTQ and yield (name, seq, qual) tuples."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline().rstrip("\n")
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # '+'
            qual = fh.readline().rstrip("\n")
            name = header[1:].split()[0]
            yield name, header, seq, qual


def orient_reads(fastq: Path, hits: dict, out_path, summary_fh=None):
    """Write oriented reads; RC those whose best adapter hit ends with '_r'."""
    counts = {"kept": 0, "rc": 0, "unstranded": 0}
    for name, header, seq, qual in open_fastq(fastq):
        target = hits.get(name, "*")
        if target == "*":
            counts["unstranded"] += 1
            strand = "."
        elif needs_rc(target):
            seq = rev_cmp(seq)
            qual = qual[::-1]
            counts["rc"] += 1
            strand = "-"
        else:
            counts["kept"] += 1
            strand = "+"
        record = f"{header}\n{seq}\n+\n{qual}\n"
        if summary_fh:
            summary_fh.write(f"{name}\t{target}\t{strand}\n")
        out_path.write(record)
    logger.info(
        "Orientation: kept=%d  rc=%d  unstranded=%d",
        counts["kept"], counts["rc"], counts["unstranded"],
    )


def main():
    args = parse_args()
    if not args.adapters and not args.tso_sequence:
        raise ValueError("Provide either --adapters or --tso_sequence.")

    with tempfile.NamedTemporaryFile(suffix=".vsearch.tsv", delete=False) as tmp:
        tsv_path = Path(tmp.name)
    with tempfile.NamedTemporaryFile(suffix=".adapters.fa", delete=False) as tmp_ad:
        tmp_adapters = Path(tmp_ad.name)

    try:
        adapters_path = args.adapters
        if adapters_path is None:
            write_adapters_fasta(args.tso_sequence, args.rt_adapter, tmp_adapters)
            adapters_path = tmp_adapters

        run_vsearch(args.fastq, adapters_path, args.min_adapter_id,
                    args.threads, tsv_path)
        hits = parse_hits(tsv_path)
        logger.info("Parsed %d read hits from vsearch.", len(hits))

        summary_fh = open(args.summary, "w") if args.summary else None
        if summary_fh:
            summary_fh.write("read_id\ttarget\tstrand\n")

        if args.output:
            # Text mode: gzip.write() rejects str when opened as "wb"
            with gzip.open(
                args.output, "wt", encoding="utf-8", newline="\n"
            ) as gz_out:
                orient_reads(args.fastq, hits, gz_out, summary_fh)
        else:
            orient_reads(args.fastq, hits, sys.stdout, summary_fh)

        if summary_fh:
            summary_fh.close()
    finally:
        tsv_path.unlink(missing_ok=True)
        tmp_adapters.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
