# @title
#!/usr/bin/env python3
"""
Strip internal node labels (e.g. "n0", "n1", ...) from a multi-tree Newick
file, keeping tip names, branch lengths, and topology intact.

Twisst calls `ete3.Tree(newick_string)` with no `format` argument, which
defaults to ete3's format=0. Under format=0, any label immediately after a
")" is interpreted as a numeric bootstrap/support value -- an alphabetic
internal node name like "n0" or "n3" is invalid and causes the whole tree
to fail to parse ("Unexpected newick format '...'"). Since twisst only
uses tree topology and tip identity (not internal support values), these
internal labels are safe to drop entirely.

Usage:
    python strip_internal_labels.py input.nwk output.nwk
"""

import re
import sys

# Matches an internal-node label: text immediately following ")" and
# preceding either ":" (has a branch length) or "," / ")" / ";" (no
# branch length, e.g. sometimes at the root).
INTERNAL_LABEL_RE = re.compile(r"\)('[^']*'|[A-Za-z0-9_.]+)([:,);])")


def strip_internal_labels(nwk: str) -> str:
    return INTERNAL_LABEL_RE.sub(lambda m: ')' + m.group(2), nwk)


def main():
    if len(sys.argv) != 3:
        print("Usage: python strip_internal_labels.py input.nwk output.nwk")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]
    n = 0
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if "\t" in line:
                og_id, nwk = line.split("\t", 1)
                fout.write(f"{og_id}\t{strip_internal_labels(nwk)}\n")
            else:
                fout.write(f"{strip_internal_labels(line)}\n")
            n += 1

    print(f"Processed {n} trees -> {out_path}")


if __name__ == "__main__":
    main()
