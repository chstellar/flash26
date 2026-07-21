#!/usr/bin/env python3
"""
Prepare and apply compactor rescue annotations for FLASH BLAST plots.

This script is intentionally small glue around the existing compactor annotation
helpers in resfungi_compactor_blast.py.  It keeps the Snakemake rules readable
while preserving the parsing and label-selection behavior that has already been
debugged for the manual rescue workflow.
"""

import argparse
import csv
import sys
from pathlib import Path


def find_repo_root(start_path):
    for parent in [start_path, *start_path.parents]:
        if (parent / "resfungi_compactor_blast.py").exists():
            return parent
    raise RuntimeError(
        "Could not locate resfungi_compactor_blast.py from "
        f"{start_path}. Run this script from the FLASH repository checkout."
    )


REPO_ROOT = find_repo_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO_ROOT))

import resfungi_compactor_blast as rescue  # noqa: E402


UNANNOTATED_LABELS = {
    "",
    "NA",
    "NAN",
    "NO MATCH",
    "NO BLAST",
    "UNANNOTATED",
    "NO PROTEIN/GENE HIT",
}


def is_unannotated_plot_row(row):
    labels = [
        row.get("Blast Label"),
        row.get("point_label"),
        row.get("direct_blast_label"),
        row.get("annotation"),
    ]
    labels = ["" if value is None else str(value).strip().upper() for value in labels]
    if "NO TARGET" in labels:
        return False
    return any(label in UNANNOTATED_LABELS for label in labels)


def read_summary_seed_rows(summary_path, anchor_len):
    rows = []
    seen = set()
    with open(summary_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for idx, row in enumerate(reader, start=1):
            if not is_unannotated_plot_row(row):
                continue
            sequence = rescue.clean_sequence_candidate(row.get("sequence"))
            if len(sequence) < anchor_len:
                continue
            anchor = sequence[-anchor_len:]
            if "N" in anchor:
                continue
            key = (
                row.get("metadata_category", ""),
                row.get("feature", ""),
                row.get("cluster", ""),
                sequence,
            )
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "seed_row": len(rows) + 1,
                    "seed_anchor": anchor,
                    "seed_extendor": sequence,
                    "metadata_category": row.get("metadata_category", ""),
                    "feature": row.get("feature", ""),
                    "cluster": row.get("cluster", ""),
                    "metadata": row.get("metadata", ""),
                    "total_samples": row.get("total_samples", ""),
                    "blast_label": row.get("Blast Label", ""),
                    "point_label": row.get("point_label", ""),
                    "summary_row": idx,
                    "seed_source": str(summary_path),
                    "raw_seed_row": "\t".join(str(row.get(col, "")) for col in fieldnames).replace(
                        "\t", "\\t"
                    ),
                }
            )
    return rows


def write_seed_outputs(rows, seeds_path, sidecar_path):
    seeds_path = Path(seeds_path)
    sidecar_path = Path(sidecar_path)
    seeds_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)

    anchors = []
    seen_anchors = set()
    for row in rows:
        anchor = row["seed_anchor"]
        if anchor not in seen_anchors:
            seen_anchors.add(anchor)
            anchors.append(anchor)

    with open(seeds_path, "w", newline="") as handle:
        for anchor in anchors:
            handle.write(f"{anchor}\n")

    columns = [
        "seed_row",
        "seed_anchor",
        "seed_extendor",
        "metadata_category",
        "feature",
        "cluster",
        "metadata",
        "total_samples",
        "blast_label",
        "point_label",
        "summary_row",
        "seed_source",
        "raw_seed_row",
    ]
    with open(sidecar_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})
    print(f"Wrote {len(anchors)} unique compactor anchors from {len(rows)} unannotated plot rows to {seeds_path}")
    print(f"Wrote compactor seed sidecar to {sidecar_path}")


def write_seed_chunks(seeds_path, chunk_outputs):
    anchors = []
    with open(seeds_path, newline="") as handle:
        anchors = [line.strip() for line in handle if line.strip()]

    chunk_paths = [Path(path) for path in chunk_outputs]
    for path in chunk_paths:
        path.parent.mkdir(parents=True, exist_ok=True)

    handles = [open(path, "w", newline="") for path in chunk_paths]
    try:
        for idx, anchor in enumerate(anchors):
            handles[idx % len(handles)].write(f"{anchor}\n")
    finally:
        for handle in handles:
            handle.close()
    print(f"Wrote {len(anchors)} anchors across {len(chunk_paths)} compactor seed chunks")


def command_merge_compactors(args):
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    with open(output, "w", newline="") as out_handle:
        writer = None
        for path_text in args.inputs:
            path = Path(path_text)
            if not path.exists() or path.stat().st_size == 0:
                continue
            with open(path, newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if not reader.fieldnames:
                    continue
                if writer is None:
                    writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=reader.fieldnames)
                    writer.writeheader()
                    wrote_header = True
                for row in reader:
                    writer.writerow({col: row.get(col, "") for col in writer.fieldnames})
        if not wrote_header:
            out_handle.write(
                "anchor\tcompactor\tid\tparent_id\tsupport\texact_support\textender_specificity\t"
                "extender_shift\ttotal_length\tnum_extended\texpected_read_count\n"
            )
    print(f"Merged {len(args.inputs)} compactor chunk TSVs into {output}")


def read_seed_sidecar(path):
    rows = []
    if not Path(path).exists():
        return rows
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader, start=1):
            row.setdefault("seed_row", idx)
            row.setdefault("seed_source", str(path))
            row.setdefault("raw_seed_row", "")
            rows.append(row)
    return rows


def command_seeds(args):
    rows = read_summary_seed_rows(args.summary, args.anchor_len)
    write_seed_outputs(rows, args.seeds, args.sidecar)
    if args.chunk_outputs:
        write_seed_chunks(args.seeds, args.chunk_outputs)


def command_select(args):
    thresholds = [float(item) for item in args.thresholds.split(",") if item.strip()]
    records = rescue.read_compactor_table(Path(args.compactors), thresholds)
    rescue.write_fasta(records, Path(args.fasta))
    rescue.write_selected_compactors(records, Path(args.selected))
    print(f"Wrote {len(records)} selected compactors to {args.selected}")
    print(f"Wrote selected compactor FASTA to {args.fasta}")


def command_annotate(args):
    selected = rescue.read_selected_compactors(Path(args.selected))
    annotations = rescue.merge_annotation_maps(
        rescue.read_annotation_table(Path(args.blastp), "blastp", "restricted"),
        rescue.read_annotation_table(Path(args.blast), "blast", "restricted"),
        rescue.read_annotation_table(Path(args.reblastp), "blastp", "unrestricted"),
        rescue.read_annotation_table(Path(args.reblast), "blast", "unrestricted"),
    )

    rescue.write_compactor_annotation_summary(selected, annotations, Path(args.compactor_annotations))
    seed_rows = read_seed_sidecar(args.sidecar)
    rescue.write_seed_annotation_summary(
        seed_rows,
        rescue.records_by_anchor(selected),
        annotations,
        Path(args.seed_annotations),
    )

    compactor_map = rescue.read_seed_compactor_annotation_map(Path(args.seed_annotations))
    anchor_map = rescue.read_compactor_anchor_annotation_map(Path(args.compactor_annotations))
    rescue.add_anchor_hits(anchor_map, rescue.build_anchor_annotation_map(compactor_map, args.anchor_len))

    rescue.fill_plot_summary_tsv(
        Path(args.plot_summary),
        Path(args.output_summary),
        compactor_map,
        anchor_map,
        args.anchor_len,
    )
    summary_key_map, summary_sequence_map, summary_rows = rescue.build_summary_compactor_maps(
        Path(args.output_summary)
    )
    rescue.fill_plot_annotation_tsv(
        Path(args.plot_blastp_annotated),
        Path(args.output_blastp_annotated),
        compactor_map,
        anchor_map,
        args.anchor_len,
        "blastp",
        summary_key_map,
        summary_sequence_map,
        summary_rows,
    )
    rescue.fill_plot_annotation_tsv(
        Path(args.plot_blast_annotated),
        Path(args.output_blast_annotated),
        compactor_map,
        anchor_map,
        args.anchor_len,
        "blast",
        summary_key_map,
        summary_sequence_map,
        summary_rows,
    )
    print(f"Wrote compactor-filled BLASTP annotations to {args.output_blastp_annotated}")
    print(f"Wrote compactor-filled BLASTN annotations to {args.output_blast_annotated}")
    print(f"Wrote compactor-filled plot summary to {args.output_summary}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    seeds = subparsers.add_parser("seeds", help="Select unannotated plot extendors and write compactor seeds.")
    seeds.add_argument("--summary", required=True)
    seeds.add_argument("--seeds", required=True)
    seeds.add_argument("--sidecar", required=True)
    seeds.add_argument("--anchor_len", type=int, default=31)
    seeds.add_argument("--chunk_outputs", nargs="*", default=[])
    seeds.set_defaults(func=command_seeds)

    merge = subparsers.add_parser("merge_compactors", help="Merge chunked compactor TSV outputs.")
    merge.add_argument("--output", required=True)
    merge.add_argument("--inputs", nargs="+", required=True)
    merge.set_defaults(func=command_merge_compactors)

    select = subparsers.add_parser("select", help="Select one representative compactor per anchor.")
    select.add_argument("--compactors", required=True)
    select.add_argument("--fasta", required=True)
    select.add_argument("--selected", required=True)
    select.add_argument("--thresholds", default="1000,100,5")
    select.set_defaults(func=command_select)

    annotate = subparsers.add_parser("annotate", help="Fill plot inputs using compactor BLAST annotations.")
    annotate.add_argument("--selected", required=True)
    annotate.add_argument("--sidecar", required=True)
    annotate.add_argument("--blast", required=True)
    annotate.add_argument("--blastp", required=True)
    annotate.add_argument("--reblast", required=True)
    annotate.add_argument("--reblastp", required=True)
    annotate.add_argument("--plot_blast_annotated", required=True)
    annotate.add_argument("--plot_blastp_annotated", required=True)
    annotate.add_argument("--plot_summary", required=True)
    annotate.add_argument("--output_blast_annotated", required=True)
    annotate.add_argument("--output_blastp_annotated", required=True)
    annotate.add_argument("--output_summary", required=True)
    annotate.add_argument("--compactor_annotations", required=True)
    annotate.add_argument("--seed_annotations", required=True)
    annotate.add_argument("--anchor_len", type=int, default=31)
    annotate.set_defaults(func=command_annotate)
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
