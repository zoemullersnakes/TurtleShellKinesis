#!/usr/bin/env python3
"""
phyloacc_summarize.py
=====================
Turns PhyloAcc's combined results (phyloacc_post.py -> results/) into:

  <out>_targets.tsv        elements passing the BF thresholds (BF1 and BF2),
                           with coordinates, BFs, expected number of
                           acceleration events in each kinetic origin, and
                           the nearest / overlapping gene
  <out>_convergent.tsv     the subset with >= MIN_ORIGINS expected independent
                           kinetic origins of acceleration  (main result)
  <out>_branch_counts.tsv  per branch: elements with P(accelerated) >= PP
                           under M2 (no branch is privileged in M2)
  <out>_null.tsv           convergence test: elements with >= MIN_ORIGINS
                           expected origins, kinetic lineages vs every
                           structurally matched set of akinetic lineages
  <out>_trio_null.tsv      three-way tip convergence test as in Sackton et al.
                           2019: elements with all three kinetic tips accelerated
                           and all other akinetic tips conserved, vs every
                           akinetic tip trio without sister pairs
  <out>_summary.txt        all counts and the parameters used

Follows Sackton et al. 2019 (ratite CNEE code, Oct 2018):
  * filter: log BF1 >= 10 and log BF2 >= 1
  * posteriors from the full model M2 (rate_postZ_M2, column suffix _3 =
    accelerated)
  * PhyloAcc's acceleration is irreversible down the tree, so the posterior
    probability that acceleration ARISES on a branch is
        P(acc on branch) - P(acc on its parent branch)
    Summed over the branches of a lineage (stem + everything below it) this
    is the expected number of acceleration events in that lineage; each
    lineage is capped at 1 ("tto"), and the sum over lineages is the
    expected number of independent origins (their floss_cl_pp).
  * candidates (their "ratite-accelerated elements"): BF1 >= 10, BF2 >= 1
    AND expected number of acceleration events on non-target branches < 1
    (their "losses in internal paleognath + tinamou branches < 1"; here every
    akinetic branch of the tree, internal and tip)
  * convergent: candidate with expected independent origins >= 1.8
A separate column counts akinetic TIPS with P(acc) >= 0.90 (their
neo_tip_loss / tin_tip_loss QC).

Kinetic lineages (independent gains of shell kinesis):
  Cuora      the Cuora trifasciata tip
  Terrapene  the Terrapene triunguis tip
  Emydinae   the Emys + Emydoidea clade (stem and both tips)

Matched null: one akinetic two-tip clade + two akinetic tips outside it,
i.e. the same shape as the kinetic set, scored the same way on all elements
(no BF filter, because BF1/BF2 are defined for the kinetic targets only).
p = (1 + #null sets with >= observed) / (1 + #null sets).

USAGE
  python phyloacc_summarize.py \
      --results /anvil/scratch/x-zlm9173/phyloacc_out/results \
      --gtf     /anvil/scratch/x-zlm9173/Chrysemys_picta_bellii.gtf \
      --tree    /anvil/scratch/x-zlm9173/phast_turtles/neutral.labeled.mod \
      --out     /anvil/scratch/x-zlm9173/phyloacc_summary/turtles
"""

import argparse
import bisect
import itertools
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# ─────────────── pre-specified parameters ───────────────
BF1_MIN = 10.0        # log BF, M1 vs M0
BF2_MIN = 1.0         # log BF, M1 vs M2
PP_MIN = 0.90         # P(accelerated) for per-branch counts / akinetic-tip QC
MIN_ORIGINS = 1.8     # convergent = expected independent origins >= this (Sackton: 1.8)
AKIN_MAX = 1.0        # candidates: expected events on akinetic branches < this

CUORA_NAMES = ["Cuotri", "Cuora_trifasciata"]
TERRAPENE = "Terrapene_triunguis"
EMYS_PAIR = ("Emydoidea_blandingii", "Emys_orbicularis")
# akinetic two-tip clades available for the matched null
NULL_PAIRS = [("Mauremys_mutica", "Mauremys_reevesii"),
              ("Graptemys_geographica", "Malaclemys_terrapin_pileata")]


def log(msg, fh):
    print(msg)
    fh.write(msg + "\n")


def read_table(path):
    """PhyloAcc tables: '#' comment lines, then a tab-separated header."""
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def find_col(cols, *patterns, required=True, what=""):
    for p in patterns:
        hits = [c for c in cols if re.fullmatch(p, c, flags=re.I)]
        if hits:
            return hits[0]
    if required:
        sys.exit(f"ERROR: no {what} column in elem_lik.txt; columns are:\n  " + "\n  ".join(cols))
    return None


def read_rates(path):
    """Locus ID + the accelerated-state (_3) column of every branch."""
    hdr = pd.read_csv(path, sep="\t", comment="#", nrows=0).columns
    keep = [hdr[0]] + [c for c in hdr if c.endswith("_3")]
    df = pd.read_csv(path, sep="\t", comment="#", usecols=keep, low_memory=False)
    df = df.rename(columns={hdr[0]: "locus"})
    df.columns = ["locus"] + [c[:-2] for c in df.columns[1:]]
    df["locus"] = df["locus"].astype(str)
    return df.set_index("locus")


def pair_label(a, b):
    return "-".join(sorted([a, b]))


class Node:
    def __init__(self):
        self.children, self.name, self.parent = [], "", None

    def tips(self):
        return [self.name] if not self.children else \
            [t for c in self.children for t in c.tips()]


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                c = node()
                c.parent = n
                n.children.append(c)
                if s[pos] == ",":
                    pos += 1
                    continue
                pos += 1
                break
        st = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        n.name = s[st:pos].strip().strip("'\"")
        if pos < len(s) and s[pos] == ":":
            pos += 1
            while pos < len(s) and s[pos] not in ",();":
                pos += 1
        return n

    return node()


def read_tree(path):
    """Newick file, or a .mod file with a TREE: line. Internal nodes are
    (re)labelled the way PhyloAcc does: sorted first tips of the two children."""
    txt = open(path).read()
    m = re.search(r"^TREE:\s*(.+;)", txt, re.M)
    root = parse_newick(m.group(1) if m else txt[txt.index("("):])

    def walk(n):
        for c in n.children:
            walk(c)
        if n.children:
            n.name = "-".join(sorted(sorted(c.tips())[0] for c in n.children))
    walk(root)
    parent, below, nodes = {}, {}, []

    def collect(n):
        nodes.append(n)
        for c in n.children:
            parent[c.name] = n.name
            collect(c)
        below[n.name] = [x.name for x in all_nodes(n)]
    collect(root)
    return root, parent, below


def all_nodes(n):
    out = [n]
    for c in n.children:
        out += all_nodes(c)
    return out


def expected_events(rates, clade_branches, parent):
    """sum over branches of P(acc) - P(acc on parent); root parent = 0."""
    tot = 0
    for b in clade_branches:
        p = parent.get(b)
        par = rates[p] if p in rates.columns else 0.0
        tot = tot + (rates[b] - par)
    return tot


def expected_origins(rates, lineages, parent):
    """lineages: list of branch lists; each capped at 1, then summed."""
    per = [expected_events(rates, br, parent).clip(upper=1) for br in lineages]
    return per, sum(per)


def read_genes(gtf, coding_only):
    genes = defaultdict(list)
    n_all = 0
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            n_all += 1
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            bt = a.get("gene_biotype", a.get("gene_type"))
            if coding_only and bt is not None and bt != "protein_coding":
                continue
            name = a.get("gene_name") or a.get("gene") or a.get("gene_id", "NA")
            genes[f[0]].append((int(f[3]) - 1, int(f[4]), name))
    for c in genes:
        genes[c].sort()
    return genes, n_all


def annotate(chrom, s, e, genes):
    """overlapping gene(s), else nearest gene up- and downstream (by coordinate)."""
    g = genes.get(chrom, [])
    if not g:
        return "NA", "no genes on sequence", "", ""
    starts = [x[0] for x in g]
    inside = [x[2] for x in g[max(0, bisect.bisect_left(starts, s) - 200):
                               bisect.bisect_left(starts, e)] if x[0] < e and x[1] > s]
    if inside:
        return ";".join(dict.fromkeys(inside)), "intronic (within gene)", 0, ""
    i = bisect.bisect_left(starts, s)
    left = [x for x in g[max(0, i - 200):i] if x[1] <= s]
    left = max(left, key=lambda x: x[1]) if left else None
    right = g[i] if i < len(g) else None
    cand = []
    if left:
        cand.append((s - left[1], left[2], "left"))
    if right:
        cand.append((right[0] - e, right[2], "right"))
    d, name, _ = min(cand)
    flank = " | ".join(f"{n} ({dd} bp {side})" for dd, n, side in sorted(cand))
    return name, "intergenic", d, flank


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="phyloacc_post results/ folder")
    ap.add_argument("--tree", required=True,
                    help="species tree (Newick, or the .mod file) used for PhyloAcc")
    ap.add_argument("--gtf", help="reference GTF (uncompressed); omit to skip genes")
    ap.add_argument("--all-genes", action="store_true",
                    help="use every gene, not only protein-coding ones")
    ap.add_argument("--bf1", type=float, default=BF1_MIN)
    ap.add_argument("--bf2", type=float, default=BF2_MIN)
    ap.add_argument("--pp", type=float, default=PP_MIN)
    ap.add_argument("--min-origins", type=float, default=MIN_ORIGINS)
    ap.add_argument("--akinetic-max", type=float, default=AKIN_MAX,
                    help="candidates need expected events on akinetic branches < this")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    lf = open(a.out + "_summary.txt", "w")
    R = a.results

    # ── element likelihoods / BFs ──
    el = read_table(os.path.join(R, "elem_lik.txt"))
    cols = list(el.columns)
    c_bf1 = find_col(cols, r"log_?bf1", r"bf1", what="BF1")
    c_bf2 = find_col(cols, r"log_?bf2", r"bf2", what="BF2")
    c_bf3 = find_col(cols, r"log_?bf3", r"bf3", required=False)
    c_best = find_col(cols, r"best\.?fit\.?model", r"best.*model", required=False)
    el = el.astype({c: str for c in cols if el[c].dtype == object})

    rates = read_rates(os.path.join(R, "rate_postZ_M2.txt"))
    rate_ids = set(rates.index[:2000])
    # which elem_lik column is the key used in rate_postZ, which holds coordinates
    c_key = max(cols, key=lambda c: el[c].astype(str).isin(rate_ids).sum())
    if el[c_key].astype(str).isin(rate_ids).sum() == 0:
        sys.exit("ERROR: no elem_lik column matches the Locus IDs of rate_postZ_M2.txt")
    coord_re = re.compile(r"^(.+):(\d+)-(\d+)$")
    c_coord = next((c for c in cols if el[c].astype(str).head(50).str.match(coord_re).all()), None)
    el["key"] = el[c_key].astype(str)
    log(f"elem_lik.txt: {len(el)} elements; key = '{c_key}', coordinates = '{c_coord}', "
        f"BF columns = {c_bf1}, {c_bf2}, {c_bf3}", lf)

    # ── tree, kinetic lineages ──
    root, parent, below = read_tree(a.tree)
    branches = list(rates.columns)
    tree_br = [n.name for n in all_nodes(root) if n is not root]
    miss = sorted(set(tree_br) - set(branches))
    if miss:
        sys.exit("ERROR: tree branches missing from rate_postZ_M2 columns: " + ", ".join(miss)
                 + "\nrate_postZ columns: " + ", ".join(branches))
    cuora = next((n for n in CUORA_NAMES if n in branches), None)
    stem = pair_label(*EMYS_PAIR)
    lineages = {"Cuora": [cuora], "Terrapene": [TERRAPENE], "Emydinae": below[stem]}
    log("Kinetic lineages (branches summed): " +
        "; ".join(f"{o} = {', '.join(b)}" for o, b in lineages.items()), lf)
    log("Parents: " + "; ".join(f"{b} <- {parent[b]}" for br in lineages.values() for b in br), lf)
    log(f"Thresholds: log BF1 >= {a.bf1}, log BF2 >= {a.bf2}, expected origins >= "
        f"{a.min_origins} (M2 posteriors, Sackton et al. 2019); P(acc) >= {a.pp} for "
        f"branch counts and akinetic-tip QC", lf)

    per, tot = expected_origins(rates, list(lineages.values()), parent)
    kin_tips = {cuora, TERRAPENE, *EMYS_PAIR}
    tips = [b for b in branches if b in set(root.tips())]
    akin_tips = [t for t in tips if t not in kin_tips]
    akin_acc = (rates[akin_tips] >= a.pp).sum(axis=1)
    kin_branches = {b for br in lineages.values() for b in br}
    akin_branches = [b for b in tree_br if b not in kin_branches]
    akin_events = expected_events(rates, akin_branches, parent)

    # ── candidate table ──
    el = el.set_index("key")
    bf1 = pd.to_numeric(el[c_bf1], errors="coerce")
    bf2 = pd.to_numeric(el[c_bf2], errors="coerce")
    tgt = el[(bf1 >= a.bf1) & (bf2 >= a.bf2)].copy()
    n_bf1 = int((bf1 >= a.bf1).sum())
    out = pd.DataFrame(index=tgt.index)
    if c_coord:
        m = tgt[c_coord].astype(str).str.extract(coord_re)
        out["chrom"], out["start"], out["end"] = m[0], m[1].astype(int), m[2].astype(int)
        out["length"] = out["end"] - out["start"]
    out["element"] = tgt[c_coord] if c_coord else tgt.index
    out["logBF1"], out["logBF2"] = bf1[tgt.index], bf2[tgt.index]
    if c_bf3:
        out["logBF3"] = pd.to_numeric(tgt[c_bf3], errors="coerce")
    if c_best:
        out["best_model"] = tgt[c_best]
    for o, v in zip(lineages, per):
        out[f"exp_events_{o}"] = v.reindex(out.index).round(3)
    out["expected_origins"] = tot.reindex(out.index).round(3)
    out["origins_with_event_ge_0.5"] = pd.concat(
        [(v.reindex(out.index) >= 0.5).rename(o) for o, v in zip(lineages, per)], axis=1
    ).apply(lambda r: ",".join(o for o in lineages if r[o]), axis=1)
    out["exp_events_akinetic"] = akin_events.reindex(out.index).round(3)
    out["akinetic_tips_accelerated"] = akin_acc.reindex(out.index)
    out["candidate"] = out["exp_events_akinetic"] < a.akinetic_max
    out["convergent"] = out["candidate"] & (out["expected_origins"] >= a.min_origins)

    if a.gtf and c_coord:
        genes, n_all = read_genes(a.gtf, not a.all_genes)
        log(f"GTF: {n_all} genes, {sum(len(v) for v in genes.values())} used "
            f"({'all' if a.all_genes else 'protein-coding where annotated'})", lf)
        ann = [annotate(c, s, e, genes) for c, s, e in zip(out["chrom"], out["start"], out["end"])]
        out["nearest_gene"], out["location"], out["distance_bp"], out["flanking_genes"] = zip(*ann) \
            if ann else ([], [], [], [])

    out = out.sort_values(["convergent", "candidate", "expected_origins", "logBF1"],
                          ascending=[False, False, False, False])
    out.to_csv(a.out + "_targets.tsv", sep="\t", index=False)
    conv = out[out["convergent"]]
    conv.to_csv(a.out + "_convergent.tsv", sep="\t", index=False)

    log(f"\nElements: {len(el)}", lf)
    log(f"  log BF1 >= {a.bf1}: {n_bf1}", lf)
    log(f"  log BF1 >= {a.bf1} and log BF2 >= {a.bf2}: {len(out)}", lf)
    log(f"  ... and expected events on akinetic branches < {a.akinetic_max} "
        f"(CANDIDATES): {int(out['candidate'].sum())}", lf)
    log("  (the counts below are for all elements passing the BF filters)", lf)
    eo = out["expected_origins"]
    log(f"    expected origins < 0.5: {(eo < 0.5).sum()};  0.5-1.5: {((eo >= 0.5) & (eo < 1.5)).sum()};"
        f"  1.5-2: {((eo >= 1.5) & (eo < 2)).sum()};  >= 2: {(eo >= 2).sum()};  >= 2.5: {(eo >= 2.5).sum()}", lf)
    for o in lineages:
        log(f"    expected events in {o} >= 0.5: {(out[f'exp_events_{o}'] >= 0.5).sum()}", lf)
    log(f"    with no akinetic tip accelerated: {(out['akinetic_tips_accelerated'] == 0).sum()}", lf)
    cand = out[out["candidate"]]
    log(f"  candidates by lineage (expected events >= 0.5): " + ", ".join(
        f"{o}: {(cand[f'exp_events_{o}'] >= 0.5).sum()}" for o in lineages), lf)
    log(f"  closest to convergent: max expected origins among candidates = "
        f"{cand['expected_origins'].max() if len(cand) else float('nan')}", lf)
    log(f"  convergent (candidate and expected origins >= {a.min_origins}): {len(conv)}  "
        f"({(conv['akinetic_tips_accelerated'] == 0).sum()} with no akinetic tip accelerated)", lf)
    if "location" in out:
        log("  convergent by location: " +
            ", ".join(f"{k}: {v}" for k, v in conv["location"].value_counts().items()), lf)

    # ── per-branch counts under M2 ──
    bc = (rates[tree_br] >= a.pp).sum().rename("elements_accelerated_M2").to_frame()
    kin_br = {b for br in lineages.values() for b in br}
    bc["kinetic"] = bc.index.isin(kin_br)
    bc = bc.sort_values("elements_accelerated_M2", ascending=False)
    bc.to_csv(a.out + "_branch_counts.tsv", sep="\t", index_label="branch")
    log(f"\nPer-branch counts (M2, P >= {a.pp}) written; kinetic branches:", lf)
    for b, r in bc[bc["kinetic"]].iterrows():
        rank = list(bc.index).index(b) + 1
        log(f"  {b}: {r['elements_accelerated_M2']}  (rank {rank} of {len(bc)})", lf)

    # ── matched null (same scoring, all elements) ──
    sets = []
    for pair in NULL_PAIRS:
        st = pair_label(*pair)
        if st not in below:
            log(f"  (null pair {pair} skipped: not a clade in the tree)", lf)
            continue
        others = [t for t in akin_tips if t not in pair]
        for t1, t2 in itertools.combinations(sorted(others), 2):
            sets.append(([below[st], [t1], [t2]], f"({pair[0]},{pair[1]}) + {t1} + {t2}"))

    def n_conv(groups):
        return int((expected_origins(rates, groups, parent)[1] >= a.min_origins).sum())

    obs = int((tot >= a.min_origins).sum())
    null = [(name, n_conv(g)) for g, name in sets]
    nd = pd.DataFrame(null, columns=["null_set", "elements_convergent_M2"])
    nd.loc[len(nd)] = ["KINETIC (observed)", obs]
    nd.to_csv(a.out + "_null.tsv", sep="\t", index=False)
    if null:
        vals = np.array([v for _, v in null])
        p = (1 + (vals >= obs).sum()) / (1 + len(vals))
        log(f"\nConvergence test (all elements, expected origins >= {a.min_origins}):", lf)
        log(f"  kinetic: {obs} elements", lf)
        log(f"  {len(vals)} matched akinetic sets: median {np.median(vals):.0f}, "
            f"range {vals.min()}-{vals.max()}", lf)
        log(f"  p = {p:.3g}  (rank of the kinetic set among all sets)", lf)
    # ── three-way tip convergence test (Sackton et al. 2019) ──
    # A trio is convergent for an element when all three tips have
    # P(acc) >= PP and every tip outside the trio's group has P(acc) < PP.
    # Kinetic trios: Cuora + Terrapene + one Emydinae tip (the other Emydinae
    # tip is left unconstrained, as other ratites were). Null trios: every
    # trio of akinetic tips with no sister pair; for them every other tip,
    # kinetic ones included, must be non-accelerated (as in their code).
    acc_tip = rates[tips] >= a.pp
    n_acc_tip = acc_tip.sum(axis=1)

    def trio_count(trio, group):
        inside = acc_tip[list(trio)].all(axis=1)
        outside = n_acc_tip - acc_tip[list(group)].sum(axis=1)
        return int((inside & (outside == 0)).sum())

    sisters = set()
    for t in tips:
        sib = [x for x in tips if x != t and parent.get(x) == parent.get(t)]
        for x in sib:
            sisters.add(frozenset((t, x)))
    kin_group = [cuora, TERRAPENE, *EMYS_PAIR]
    kin_trios = [(cuora, TERRAPENE, e) for e in EMYS_PAIR]
    obs_trio = [(", ".join(tr), trio_count(tr, kin_group)) for tr in kin_trios]
    null_trios = [tr for tr in itertools.combinations(sorted(akin_tips), 3)
                  if not any(frozenset(pr) in sisters for pr in itertools.combinations(tr, 2))]
    null_trio = [(", ".join(tr), trio_count(tr, tr)) for tr in null_trios]
    td = pd.DataFrame(null_trio, columns=["trio", "elements_3way_convergent"])
    td["kinetic"] = False
    td = pd.concat([pd.DataFrame([(n, v, True) for n, v in obs_trio], columns=td.columns), td])
    td.to_csv(a.out + "_trio_null.tsv", sep="\t", index=False)
    if null_trio:
        nv = np.array([v for _, v in null_trio])
        log(f"\nThree-way tip convergence (Sackton et al. 2019; all elements, M2, P >= {a.pp}; "
            f"all three tips accelerated, all other non-target tips conserved):", lf)
        for name, v in obs_trio:
            p_t = (1 + (nv >= v).sum()) / (1 + len(nv))
            log(f"  kinetic trio {name}: {v} elements  (p = {p_t:.3g})", lf)
        log(f"  {len(nv)} akinetic trios (no sister pairs): median {np.median(nv):.0f}, "
            f"mean {nv.mean():.1f}, range {nv.min()}-{nv.max()}", lf)
        prop = [v / max(1, int((acc_tip[list(tr)].all(axis=1)).sum())) for tr in kin_trios]
        log(f"  (elements with all three kinetic tips accelerated, ignoring the rest of the tree: "
            + ", ".join(str(int(acc_tip[list(tr)].all(axis=1).sum())) for tr in kin_trios) + ")", lf)
    log(f"\nWrote {a.out}_targets.tsv, _convergent.tsv, _branch_counts.tsv, _null.tsv, "
        f"_trio_null.tsv, _summary.txt", lf)


if __name__ == "__main__":
    main()