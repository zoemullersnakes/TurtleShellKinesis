#!/usr/bin/env python3
"""
Move (or copy) the alignment FASTA files for a list of candidate orthogroups
into their own folder, based on a twisst/RERconverge TSV that has a 'name'
column of OG IDs.

Edit the four variables below, then run:
    python3 move_twisst_candidates.py
"""

import shutil
from pathlib import Path
import pandas as pd

# ---- EDIT THESE ----------------------------------------------------------
TSV_PATH        = "twisst_putative_genes.tsv"        # your candidate list
ALIGNMENT_DIR   = "/Users/zlm9173/Desktop/Turtles/alignments"   # <-- where ALL OG alignments currently live (fix this path)
DEST_DIR        = "/Users/zlm9173/Desktop/Turtles/twisst_candidates"  # where to put the matched ones
FILENAME_PATTERN = "{og}_aligned.fasta"               # matches OG0012875_aligned.fasta style naming
ACTION           = "copy"                             # "copy" (safe, default) or "move"
# ---------------------------------------------------------------------------

def main():
    df = pd.read_csv(TSV_PATH, sep=",")
    if "name" not in df.columns:
        raise SystemExit(f"No 'name' column found in {TSV_PATH}. Columns are: {list(df.columns)}")

    og_ids = df["name"].dropna().unique().tolist()
    print(f"Found {len(og_ids)} orthogroup IDs in {TSV_PATH}")

    src_dir = Path(ALIGNMENT_DIR)
    dest_dir = Path(DEST_DIR)
    dest_dir.mkdir(parents=True, exist_ok=True)

    found, missing = [], []

    for og in og_ids:
        src_file = src_dir / FILENAME_PATTERN.format(og=og)
        if src_file.exists():
            dest_file = dest_dir / src_file.name
            if ACTION == "move":
                shutil.move(str(src_file), str(dest_file))
            else:
                shutil.copy2(str(src_file), str(dest_file))
            found.append(og)
        else:
            missing.append(og)

    print(f"\n{ACTION.capitalize()}d {len(found)} alignment files to {dest_dir}")
    if missing:
        print(f"\n{len(missing)} orthogroups had no matching file in {src_dir} "
              f"(check FILENAME_PATTERN / ALIGNMENT_DIR):")
        for og in missing:
            print(f"  - {FILENAME_PATTERN.format(og=og)}")

if __name__ == "__main__":
    main()