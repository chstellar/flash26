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


def first_present(row, names):
    for name in names:
        value = row.get(name)
        if value is not None and value.strip() and value.strip().upper() != "NA":
            return value.strip()
    return ""


def seed_anchor(seed, anchor_len):
    if anchor_len <= 0:
        raise ValueError("--anchor_len must be positive")
    return seed[-anchor_len:]


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
    parser.add_argument(
        "--target_len",
        type=int,
        default=None,
        help="Optional target length, used only for sanity reporting.",
    )
    args = parser.parse_args()

    exact = {}
    for row in read_tsv(Path(args.seed_annotations)):
        seed = first_present(row, ["seed_extendor", "extendor", "seed", "sequence", "query"])
        compactor = first_present(row, ["compactor_sequence", "compactor"])
        if seed and compactor:
            exact.setdefault(seed, compactor)

    by_anchor = {}
    if args.selected:
        for row in read_tsv(Path(args.selected)):
            anchor = first_present(row, ["anchor"])
            compactor = first_present(row, ["compactor_sequence", "compactor"])
            if anchor and compactor:
                by_anchor.setdefault(anchor, compactor)

    missing = 0
    short = 0
    expected_seed_len = None
    if args.target_len is not None:
        expected_seed_len = args.anchor_len + args.target_len

    with open(args.output, "w", newline="") as out:
        for seed in read_lines(Path(args.input)):
            if expected_seed_len is not None and len(seed) < expected_seed_len:
                short += 1
            compactor = exact.get(seed)
            if compactor is None and args.selected:
                compactor = by_anchor.get(seed_anchor(seed, args.anchor_len))
            if compactor is None:
                missing += 1
                compactor = "NA"
            out.write(compactor + "\n")

    print(f"Wrote {args.output}; missing compactors for {missing} seed(s).")
    if short:
        print(
            f"Warning: {short} input seed(s) were shorter than "
            f"anchor_len + target_len = {expected_seed_len}."
        )


if __name__ == "__main__":
    main()
