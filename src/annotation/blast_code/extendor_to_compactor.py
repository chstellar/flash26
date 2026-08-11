#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_lines(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main():
    parser = argparse.ArgumentParser(
        description="Map one extendor/seed sequence per line to its compactor sequence."
    )
    parser.add_argument("--input", default="test.txt")
    parser.add_argument("--output", default="compactor.txt")
    parser.add_argument(
        "--seed_annotations",
        required=True,
        help="*_compactor_fungus_regular_seed_annotations.tsv from compactor_reannotation.",
    )
    parser.add_argument(
        "--selected",
        default=None,
        help="Optional *_compactor_fungus_regular_selected.tsv fallback, matched by anchor suffix.",
    )
    parser.add_argument("--anchor_len", type=int, default=31)
    args = parser.parse_args()

    exact = {}
    for row in read_tsv(Path(args.seed_annotations)):
        seed = (row.get("seed_extendor") or "").strip()
        compactor = (
            row.get("compactor_sequence")
            or row.get("compactor")
            or ""
        ).strip()
        if seed and compactor and compactor.upper() != "NA":
            exact.setdefault(seed, compactor)

    by_anchor = {}
    if args.selected:
        for row in read_tsv(Path(args.selected)):
            anchor = (row.get("anchor") or "").strip()
            compactor = (
                row.get("compactor_sequence")
                or row.get("compactor")
                or ""
            ).strip()
            if anchor and compactor and compactor.upper() != "NA":
                by_anchor.setdefault(anchor, compactor)

    missing = 0
    with open(args.output, "w", newline="") as out:
        for seed in read_lines(Path(args.input)):
            compactor = exact.get(seed)
            if compactor is None and args.selected:
                compactor = by_anchor.get(seed[-args.anchor_len:])
            if compactor is None:
                missing += 1
                compactor = "NA"
            out.write(compactor + "\n")

    print(f"Wrote {args.output}; missing compactors for {missing} seed(s).")


if __name__ == "__main__":
    main()
