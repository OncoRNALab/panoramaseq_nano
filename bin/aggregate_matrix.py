#!/usr/bin/env python3
"""Aggregate HDF matrix chunks from create_matrix and write 10x MEX output."""
import argparse
import logging
import sys
from pathlib import Path

_bin_dir = Path(__file__).resolve().parent
if str(_bin_dir) not in sys.path:
    sys.path.insert(0, str(_bin_dir))

from workflow_glue.expression_matrix import ExpressionMatrix  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [aggregate_matrix] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("aggregate_matrix")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("hdf_dir", type=Path, help="Directory with *.gene.hdf / *.transcript.hdf chunks.")
    p.add_argument("output_dir", type=Path, help="Parent output directory for MEX folders.")
    p.add_argument("--prefix", type=str, default="sample", help="Sample prefix for folder names.")
    p.add_argument(
        "--features", nargs="+", default=["gene", "transcript"],
        choices=["gene", "transcript"],
        help="Feature types to aggregate (default: gene transcript).",
    )
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for feature in args.features:
        hdfs = sorted(args.hdf_dir.glob(f"*.{feature}.hdf"))
        if not hdfs:
            logger.warning("No %s HDF chunks found in %s; skipping.", feature, args.hdf_dir)
            continue
        logger.info("Aggregating %d %s HDF chunk(s).", len(hdfs), feature)
        matrix = ExpressionMatrix.aggregate_hdfs(hdfs)
        out = args.output_dir / f"{args.prefix}_{feature}_bc_matrix"
        matrix.to_mex(str(out))
        logger.info("Wrote %s", out)


if __name__ == "__main__":
    main()
