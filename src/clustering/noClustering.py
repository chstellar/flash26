# noClustering.py
# Daniel Cotter
# this script is used to bypass clustering for short anchors
# it takes as input a tsv with anchors and outputs a tsv with a
# line number per anchor
# column 1 is the line number, column 2 is the anchor
# is also takes a temp_dir paramater that is unused but
# is required by other scripts in the snakemake rule

import pandas as pd
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Process anchors without clustering.")
    parser.add_argument(
        "--input",
        type=str,
        help="Input file containing anchors in tsv format.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file to write anchors with line numbers.",
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        help="Temporary directory for output files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    temp_dir = args.temp_dir
    # Ensure the temp_dir exists
    Path(temp_dir).mkdir(parents=True, exist_ok=True)
    # Read the input file
    anchors_df = pd.read_csv(input_file, sep="\t", header=None, names=["anchor"])
    # Add a line number column
    anchors_df["line_number"] = range(1, len(anchors_df) + 1)
    # Reorder the columns to have line_number first
    anchors_df = anchors_df[["line_number", "anchor"]]
    # Write the output file
    anchors_df.to_csv(output_file, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
