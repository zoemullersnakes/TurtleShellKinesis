#!/usr/bin/env python3
"""
make_phyloacc_input.py
======================
Builds PhyloAcc's concatenated-alignment input (-a / -b) for conserved
NON-CODING elements, straight from the per-window MAFs of the phastCons
pipeline (already cleaned: no "# hal" comment lines, tree tips only).

Steps
  1. Read conserved elements (phastCons most-conserved BED).
  2. Drop elements that overlap any exon (--exclude BED) and elements shorter
     than --min-len. Both filters are applied as given; nothing is tuned.
  3. For every remaining element, take the alignment columns at each
     reference position (reference-projected: columns where the reference has
     a gap are not included). Species absent from a block, and reference
     positions not covered by any block, become "-".
  4. Write
       <out>.fa      one sequence per species, all elements concatenated;
                     headers = tree tips exactly; original case kept so
                     PhyloAcc --softmask can mask repeats if wanted
       <out>.bed     chrom, start, end (0-based, end-exclusive, in the
                     concatenated alignment), element ID  -> PhyloAcc -b
       <out>.elements.tsv  ID -> genomic coordinates, length, % gaps per element
       <out>.log     counts at every filtering step

USAGE
  python make_phyloacc_input.py \
      --elements phast_turtles/final/most-cons.all.bed \
      --exclude  exons.bed \
      --windows  phast_turtles/windows.txt \
      --maf-dir  phast_turtles/maf_windows \
      --tree     phast_turtles/turtles.noanc.nh \
      --ref      Chrysemys_picta_bellii \
      --min-len  50 \
      --out      phyloacc_input/cne
"""

import argparse
import bisect
import os
import re
import sys
from collections import defaultdict

import numpy as np


def log(msg, fh):
    print(msg)
    fh.write(msg + "\n")
    fh.flush()


def read_bed(path):
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            out[f[0]].append((int(f[1]), int(f[2])))
    for c in out:
        out[c].sort()
    return out


def merge(ivs):
    m = []
    for s, e in ivs:
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return m


def overlaps_any(s, e, starts, ends):
    """ivs merged & sorted: does [s,e) overlap any?"""
    i = bisect.bisect_right(starts, s) - 1
    if i >= 0 and ends[i] > s:
        return True
    return i + 1 < len(starts) and starts[i + 1] < e


def tree_tips(path):
    t = open(path).read()
    t = re.sub(r":[0-9.eE+-]+", "", t)
    return [x for x in re.split(r"[(),;\s]+", t) if x]


def maf_blocks(path):
    """yield list of (src, start, size, strand, srcsize, text) per block"""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("a"):
                if rows:
                    yield rows
                rows = []
            elif line.startswith("s"):
                f = line.split()
                rows.append((f[1], int(f[2]), int(f[3]), f[4], int(f[5]), f[6]))
            elif not line.strip() and rows:
                yield rows
                rows = []
    if rows:
        yield rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--exclude", help="BED of exons (any overlap -> element dropped)")
    ap.add_argument("--windows", required=True, help="windows.txt: ID CHR START LEN")
    ap.add_argument("--maf-dir", required=True)
    ap.add_argument("--tree", required=True, help="Newick with the species as tips")
    ap.add_argument("--ref", required=True, help="reference genome name in the MAF")
    ap.add_argument("--min-len", type=int, default=50)
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    lf = open(a.out + ".log", "w")
    species = tree_tips(a.tree)
    sp_index = {s: i for i, s in enumerate(species)}
    ref_i = sp_index.get(a.ref)
    if ref_i is None:
        sys.exit(f"ERROR: reference {a.ref} is not a tip of {a.tree}: {species}")
    log(f"Species ({len(species)}): {', '.join(species)}", lf)

    elems = read_bed(a.elements)
    n_all = sum(len(v) for v in elems.values())
    excl = {c: merge(v) for c, v in read_bed(a.exclude).items()} if a.exclude else {}
    kept, n_exon, n_short = {}, 0, 0
    for c, ivs in elems.items():
        ex = excl.get(c, [])
        xs, xe = [x[0] for x in ex], [x[1] for x in ex]
        keep = []
        for s, e in ivs:
            if ex and overlaps_any(s, e, xs, xe):
                n_exon += 1
            elif e - s < a.min_len:
                n_short += 1
            else:
                keep.append((s, e))
        if keep:
            kept[c] = keep
    n_keep = sum(len(v) for v in kept.values())
    log(f"Elements read: {n_all}", lf)
    log(f"  dropped, overlap exon:       {n_exon}", lf)
    log(f"  dropped, shorter than {a.min_len} bp: {n_short}", lf)
    log(f"  kept (non-coding, >= {a.min_len} bp): {n_keep}  "
        f"({sum(e - s for v in kept.values() for s, e in v)} bp)", lf)

    win = defaultdict(list)
    with open(a.windows) as fh:
        for line in fh:
            wid, c, s, n = line.split()[:4]
            win[c].append((int(s), wid))

    tmp = [open(f"{a.out}.{i}.tmp", "w") for i in range(len(species))]
    bed = open(a.out + ".bed", "w")
    tsv = open(a.out + ".elements.tsv", "w")
    tsv.write("id\tchrom\tstart\tend\tlength\taln_start\taln_end\tpct_gap\tn_species_present\n")
    offset = 0
    dup_rows = minus_ref = 0

    for c in sorted(kept):
        ivs = kept[c]
        starts = [s for s, _ in ivs]
        # per element: (n_species x len) array of '-'
        arrs = [np.full((len(species), e - s), ord("-"), dtype=np.uint8) for s, e in ivs]
        for _, wid in sorted(win.get(c, [])):
            path = os.path.join(a.maf_dir, wid + ".maf")
            if not os.path.exists(path):
                log(f"WARNING: missing MAF {path}", lf)
                continue
            for rows in maf_blocks(path):
                ref = None
                for r in rows:
                    g, _, seqname = r[0].partition(".")
                    if g == a.ref and seqname == c:
                        ref = r
                        break
                if ref is None:
                    continue
                if ref[3] != "+":
                    minus_ref += 1
                    continue
                bs, be = ref[1], ref[1] + ref[2]
                # elements overlapping [bs, be)
                j = bisect.bisect_right(starts, bs) - 1
                j = max(j, 0)
                hits = []
                while j < len(ivs) and ivs[j][0] < be:
                    if ivs[j][1] > bs:
                        hits.append(j)
                    j += 1
                if not hits:
                    continue
                reftext = np.frombuffer(ref[5].encode(), dtype=np.uint8)
                refcols = np.nonzero(reftext != ord("-"))[0]   # column of each ref base
                seen = set()
                texts = {}
                for r in rows:
                    g = r[0].partition(".")[0]
                    i = sp_index.get(g)
                    if i is None:
                        continue
                    if i in seen:
                        dup_rows += 1
                        continue
                    seen.add(i)
                    texts[i] = np.frombuffer(r[5].encode(), dtype=np.uint8)
                for j in hits:
                    es, ee = ivs[j]
                    os_, oe = max(bs, es), min(be, ee)
                    cols = refcols[os_ - bs:oe - bs]
                    for i, t in texts.items():
                        arrs[j][i, os_ - es:oe - es] = t[cols]
        for (s, e), arr in zip(ivs, arrs):
            L = e - s
            eid = f"{c}:{s}-{e}"
            for i in range(len(species)):
                tmp[i].write(arr[i].tobytes().decode())
            bed.write(f"{c}\t{offset}\t{offset + L}\t{eid}\n")
            gaps = (arr == ord("-"))
            present = int((~gaps).any(axis=1).sum())
            tsv.write(f"{eid}\t{c}\t{s}\t{e}\t{L}\t{offset}\t{offset + L}\t"
                      f"{100 * gaps.mean():.1f}\t{present}\n")
            offset += L
        log(f"  {c}: {len(ivs)} elements", lf)

    for f in tmp:
        f.close()
    bed.close()
    tsv.close()
    with open(a.out + ".fa", "w") as fa:
        for i, sp in enumerate(species):
            fa.write(f">{sp}\n")
            with open(f"{a.out}.{i}.tmp") as t:
                while True:
                    chunk = t.read(1 << 24)
                    if not chunk:
                        break
                    for k in range(0, len(chunk), 100):
                        fa.write(chunk[k:k + 100] + "\n")
            os.remove(f"{a.out}.{i}.tmp")
    if dup_rows:
        log(f"WARNING: {dup_rows} duplicate species rows ignored (first row kept)", lf)
    if minus_ref:
        log(f"WARNING: {minus_ref} blocks with minus-strand reference skipped", lf)
    log(f"Wrote {a.out}.fa ({offset} columns x {len(species)} species), "
        f"{a.out}.bed, {a.out}.elements.tsv", lf)


if __name__ == "__main__":
    main()