#!/usr/bin/env python3
"""
Join a twisst weights.csv file to gene/orthogroup names.

twisst's output has no concept of a "tree name" -- its rows are just in
the same order as the input tree file's lines. This script re-attaches
names using a simple line_number -> name lookup table (the kind produced
alongside a bare, twisst-ready tree file -- see strip_internal_labels.py /
the line_to_OG.tsv file from earlier in this project).

For each row it also reports:
  - the raw weight/count for every topology
  - each topology as a fraction of that row's total samples
  - which topology "wins" (highest count) for that gene, and by how much

Usage:
    python join_twisst_names.py \
        --weights output.twisst.weights.csv \
        --line_to_name mega_tree.twisst.line_to_OG.tsv \
        --out twisst_weights_named.tsv

    # Optionally restrict to a comparison between two GROUPS of topologies,
    # e.g. topologies representing a convergence hypothesis vs. the
    # topologies representing the expected species-tree relationship:
    python join_twisst_names.py \
        --weights output.twisst.weights.csv \
        --line_to_name mega_tree.twisst.line_to_OG.tsv \
        --out twisst_weights_named.tsv \
        --group_a topo1,topo2,topo3 \
        --group_b topo4,topo5,topo6 \
        --group_a_label convergence \
        --group_b_label species_tree
"""

import argparse
import csv


def load_line_to_name(path):
    mapping = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            if not row:
                continue
            line_number, name = row[0], row[1]
            mapping[int(line_number)] = name
    return mapping


def load_weights(path):
    """Returns (topo_names, data_rows) where data_rows is a list of
    lists of ints, one per tree, in file order."""
    with open(path) as f:
        lines = f.readlines()

    # Skip any leading comment lines (topology definitions, start with "#")
    start = 0
    while lines[start].startswith("#"):
        start += 1

    topo_names = lines[start].strip().split("\t")
    data_rows = []
    for line in lines[start + 1:]:
        line = line.strip()
        if not line:
            continue
        data_rows.append([int(x) for x in line.split("\t")])

    return topo_names, data_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weights", required=True, help="twisst weights.csv file")
    ap.add_argument("--line_to_name", required=True,
                     help="TSV with header line_number<TAB>name, mapping "
                          "each input tree's line number to a gene/OG name")
    ap.add_argument("--out", required=True, help="output TSV path")
    ap.add_argument("--group_a", default=None,
                     help="comma-separated topology names to sum as 'group A' "
                          "(e.g. a convergence hypothesis), optional")
    ap.add_argument("--group_b", default=None,
                     help="comma-separated topology names to sum as 'group B' "
                          "(e.g. the species-tree expectation), optional")
    ap.add_argument("--group_a_label", default="group_a")
    ap.add_argument("--group_b_label", default="group_b")
    args = ap.parse_args()

    name_map = load_line_to_name(args.line_to_name)
    topo_names, data_rows = load_weights(args.weights)

    group_a_idx = [topo_names.index(t) for t in args.group_a.split(",")] if args.group_a else None
    group_b_idx = [topo_names.index(t) for t in args.group_b.split(",")] if args.group_b else None

    with open(args.out, "w") as out:
        header = ["name"] + topo_names + ["total_samples", "top_topology", "top_topology_frac"]
        if group_a_idx:
            header += [f"{args.group_a_label}_count", f"{args.group_a_label}_frac"]
        if group_b_idx:
            header += [f"{args.group_b_label}_count", f"{args.group_b_label}_frac"]
        out.write("\t".join(header) + "\n")

        for i, row in enumerate(data_rows, start=1):
            name = name_map.get(i, f"line{i}")
            total = sum(row)
            top_i = max(range(len(row)), key=lambda k: row[k])
            top_name = topo_names[top_i]
            top_frac = row[top_i] / total if total else 0.0

            fields = [name] + [str(v) for v in row] + [str(total), top_name, f"{top_frac:.4f}"]

            if group_a_idx:
                a_count = sum(row[k] for k in group_a_idx)
                a_frac = a_count / total if total else 0.0
                fields += [str(a_count), f"{a_frac:.4f}"]
            if group_b_idx:
                b_count = sum(row[k] for k in group_b_idx)
                b_frac = b_count / total if total else 0.0
                fields += [str(b_count), f"{b_frac:.4f}"]

            out.write("\t".join(fields) + "\n")

    print(f"Wrote {len(data_rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
