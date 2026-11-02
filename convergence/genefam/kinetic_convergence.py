# @title
#!/usr/bin/env python3

import argparse
import itertools
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

# ─────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────
# CAFE5 output folder(s). In Colab/Jupyter, list every replicate here
# (one folder per run); on the command line, pass them as arguments instead.
DEFAULT_DIRS = [
    "/content/drive/MyDrive/turtles/cafe5_results_convergence",
    "/content/drive/MyDrive/turtles/cafe5_results_convergence2",
    "/content/drive/MyDrive/turtles/cafe5_results_convergence3",
    "/content/drive/MyDrive/turtles/cafe5_results_convergence4",
    "/content/drive/MyDrive/turtles/cafe5_results_convergence5"
]
MODEL    = "Base"          # file prefix: "Gamma" (k=3 run) or "Base"
# Where result tables are written. In Colab, keep this on Drive so the
# files survive the session.
OUT_DIR  = Path("/content/drive/MyDrive/turtles/kinetic_convergence_out")
P_THRESH = 0.05             # family-wide p-value
BP_THRESH = 0.05            # branch-probability screen only (exploratory)

# Replicates count as converged if -lnL differs by less than LNL_TOL
# log-likelihood units and lambda by less than LAMBDA_RTOL (relative).
LNL_TOL     = 1.0
LAMBDA_RTOL = 0.05

# Tip names exactly as in the CAFE5 input tree (without the <N> suffix).
# Several aliases may be listed per tip; the first one present is used.
KINETIC_ORIGINS = {
    "Cuora":          [["Cuotri", "Cuora_trifasciata"]],
    "Terrapene":      [["Terrapene_triunguis"]],
    "Emydoidea+Emys": [["Emydoidea_blandingii"], ["Emys_orbicularis"]],
}

# Closest AKINETIC relatives of each origin, for the sister-pair contrast.
# Actinemys is the akinetic sister of (Emydoidea, Emys); Terrapene's sister
# clade contains Actinemys, and Clemmys is the next akinetic outgroup.
SISTERS = {
    "Cuora":          ["Mauremys_reevesii", "Mauremys_mutica"],
    "Terrapene":      ["Actinemys_pallida", "Clemmys_guttata"],
    "Emydoidea+Emys": ["Actinemys_pallida"],
}

# Families with >20% optimizer failures (from *_results.txt), flagged in output.
UNSTABLE_FAMILIES = set()   # e.g. {"OG0000024", "OG0000030"}

# In Jupyter/Colab, sys.argv holds the kernel's own flags (e.g. "-f
# kernel.json"), so command-line parsing is skipped there and the
# DEFAULT_DIRS / MODEL / OUT_DIR settings above are used instead.
IN_NOTEBOOK = "ipykernel" in sys.modules
if IN_NOTEBOOK:
    RUN_DIRS = [Path(d) for d in DEFAULT_DIRS]
else:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="*", default=DEFAULT_DIRS,
                    help="one or more CAFE5 output directories (replicates)")
    ap.add_argument("--model", default=MODEL)
    args = ap.parse_args()
    MODEL = args.model
    RUN_DIRS = [Path(d) for d in args.dirs]
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# FILE LOADING
# ─────────────────────────────────────────────────────────────────
def cafe_file(run_dir, suffix, required=True):
    path = run_dir / f"{MODEL}_{suffix}"
    if not path.exists():
        if not required:
            return None
        avail = sorted(p.name for p in run_dir.glob("*_change.tab"))
        sys.exit(f"ERROR: {path} not found. Change files present in "
                 f"{run_dir}: {avail or 'none'}. Check the directory and "
                 "--model.")
    return path


def strip_label(label):
    """'Emys_orbicularis<8>*_5' -> 'Emys_orbicularis<8>'; '<16>_3' -> '<16>'."""
    m = re.match(r"^\s*([^<]*<\d+>)", label.replace("*", ""))
    return m.group(1) if m else label.strip()


def parse_asr_tree(path):
    """Read the first tree in CAFE5's *_asr.tre and return
    (parent, children, root) keyed by CAFE5 node labels such as 'Name<3>'."""
    newick = None
    for line in open(path):
        if re.match(r"^\s*TREE\s", line, re.I) and "=" in line:
            newick = line.split("=", 1)[1].strip()
            break
    if newick is None:
        sys.exit(f"ERROR: no TREE line found in {path}")

    parent, children = {}, {}
    pos = 0

    def node():
        nonlocal pos
        kids = []
        if newick[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                if newick[pos] == ",":
                    pos += 1
                    continue
                pos += 1          # ')'
                break
        start = pos
        while newick[pos] not in ",():;":
            pos += 1
        label = strip_label(newick[start:pos])
        if newick[pos] == ":":
            pos += 1
            while newick[pos] not in ",();":
                pos += 1
        children[label] = kids
        for k in kids:
            parent[k] = label
        return label

    root = node()
    return parent, children, root


def parse_results(run_dir):
    """-lnL, lambda(s) and alpha from *_results.txt (NaN if absent)."""
    out = {"lnL": float("nan"), "lambda": "", "alpha": ""}
    path = cafe_file(run_dir, "results.txt", required=False)
    if path is None:
        return out
    txt = path.read_text()
    m = re.search(r"-lnL\)?:?\s*([-\d.eE+]+)", txt)
    if m:
        out["lnL"] = float(m.group(1))
    m = re.search(r"^\s*Lambda:\s*(.+)$", txt, re.M)
    if m:
        out["lambda"] = m.group(1).strip()
    m = re.search(r"^\s*Alpha:\s*(.+)$", txt, re.M)
    if m:
        out["alpha"] = m.group(1).strip()
    return out


def load_run(run_dir):
    change = pd.read_csv(cafe_file(run_dir, "change.tab"), sep="\t", index_col=0)
    count  = pd.read_csv(cafe_file(run_dir, "count.tab"),  sep="\t", index_col=0)
    change.index = change.index.astype(str)
    count.index  = count.index.astype(str)
    fam = pd.read_csv(cafe_file(run_dir, "family_results.txt"), sep="\t")
    fam = fam.iloc[:, :3]
    fam.columns = ["FamilyID", "pvalue", "Significant"]
    fam["FamilyID"] = fam["FamilyID"].astype(str)
    fam = fam.set_index("FamilyID")
    if set(fam.index) - set(change.index):
        sys.exit(f"ERROR ({run_dir}): family_results IDs missing from "
                 "change.tab — files from different CAFE5 runs?")
    bp = None
    bp_path = cafe_file(run_dir, "branch_probabilities.tab", required=False)
    if bp_path is not None:
        bp = pd.read_csv(bp_path, sep="\t", index_col=0,
                         na_values=["N/A", "NA", "nan", "-"])
        bp.index = bp.index.astype(str).str.lstrip("#")
        bp.columns = [str(c).strip() for c in bp.columns]
        bp = bp.apply(pd.to_numeric, errors="coerce")
        bp = bp.reindex(index=change.index, columns=change.columns)
    else:
        print(f"  NOTE ({run_dir}): no {MODEL}_branch_probabilities.tab; "
              "branch-probability screen skipped for this run.")
    return {"dir": run_dir, "change": change, "count": count, "fam": fam,
            "bp": bp, **parse_results(run_dir)}


print(f"Loading {len(RUN_DIRS)} CAFE5 '{MODEL}' run(s)...")
runs = [load_run(d) for d in RUN_DIRS]

# ─────────────────────────────────────────────────────────────────
# REPLICATE CONVERGENCE CHECK
# ─────────────────────────────────────────────────────────────────
def first_float(s):
    m = re.search(r"[-\d.]+(?:[eE][-+]?\d+)?", str(s))
    return float(m.group(0)) if m else float("nan")


rep_table = pd.DataFrame([{
    "run": r["dir"].name or str(r["dir"]),
    "-lnL": r["lnL"], "lambda": r["lambda"], "alpha": r["alpha"],
    "n_sig_families": int((r["fam"]["pvalue"] < P_THRESH).sum()),
} for r in runs])

if len(runs) > 1:
    print("\n" + "=" * 70)
    print("  REPLICATE CONVERGENCE")
    print("=" * 70)
    print(rep_table.to_string(index=False))
    lnl = rep_table["-lnL"]
    lam = rep_table["lambda"].apply(first_float)
    if lnl.isna().any():
        print("\n  NOTE: -lnL not found for some runs (no *_results.txt?); "
              "using the first directory as primary.")
        best = 0
    else:
        best = int(lnl.idxmin())
        spread = lnl.max() - lnl.min()
        lam_spread = (lam.max() - lam.min()) / lam.median() if lam.notna().all() else float("nan")
        ok = spread < LNL_TOL and (pd.isna(lam_spread) or lam_spread < LAMBDA_RTOL)
        print(f"\n  -lnL spread: {spread:.3f}   lambda relative spread: "
              f"{lam_spread:.3%}")
        print("  Replicates CONVERGED." if ok else
              "  WARNING: replicates did NOT converge to the same optimum. "
              "Rerun (more replicates or different starting lambda) before "
              "trusting any single run; results below use the best run.")
    print(f"  Primary run (lowest -lnL): {rep_table.loc[best, 'run']}")
    rep_table.to_csv(OUT_DIR / "kinetic_replicates.tsv", sep="\t", index=False)
else:
    best = 0

primary = runs[best]
runs = [primary] + [r for i, r in enumerate(runs) if i != best]

# ─────────────────────────────────────────────────────────────────
# TREE (from the primary run; all replicates must share it)
# ─────────────────────────────────────────────────────────────────
parent, children, root = parse_asr_tree(cafe_file(primary["dir"], "asr.tre"))
tips = [n for n, k in children.items() if not k]

for r in runs:
    missing_cols = [n for n in parent if n not in r["change"].columns]
    if missing_cols:
        sys.exit(f"ERROR ({r['dir']}): tree nodes not found in change.tab: "
                 f"{missing_cols[:5]}. Were all replicates run on the same "
                 "tree and input?")

name_to_tip = {t.split("<")[0]: t for t in tips}


def resolve(aliases):
    for a in aliases:
        if a in name_to_tip:
            return name_to_tip[a]
    sys.exit(f"ERROR: none of {aliases} found among tree tips: "
             f"{sorted(name_to_tip)}")


def descendants_tips(n):
    return [n] if not children[n] else \
        [t for c in children[n] for t in descendants_tips(c)]


def ancestors(n):
    out = []
    while n in parent:
        out.append(n)
        n = parent[n]
    return out            # n itself up to (not incl.) root


def mrca(nodes):
    common = set(ancestors(nodes[0]))
    for n in nodes[1:]:
        common &= set(ancestors(n))
    return min(common, key=lambda x: len(descendants_tips(x)))


def make_unit(unit_tips):
    stem = mrca(unit_tips) if len(unit_tips) > 1 else unit_tips[0]
    clade = descendants_tips(stem)
    # branches from each tip up to and including the origin's stem branch
    paths = {}
    for t in unit_tips:
        up = ancestors(t)
        paths[t] = up[:up.index(stem) + 1]
    branches = sorted({b for p in paths.values() for b in p})
    return {"stem": stem, "tips": unit_tips, "clade": clade,
            "paths": paths, "branches": branches}


origins = {}
for name, tip_aliases in KINETIC_ORIGINS.items():
    u = make_unit([resolve(a) for a in tip_aliases])
    if set(u["clade"]) != set(u["tips"]):
        print(f"  WARNING: origin '{name}' MRCA also contains "
              f"{sorted(set(u['clade']) - set(u['tips']))}")
    origins[name] = u

K = len(origins)
kinetic_tips     = {t for u in origins.values() for t in u["tips"]}
kinetic_branches = {b for u in origins.values() for b in u["branches"]}
bg_tips     = [t for t in tips if t not in kinetic_tips]
bg_internal = [n for n in parent if children[n] and n not in kinetic_branches]
sister_cols = {name: [resolve([s]) for s in sis] for name, sis in SISTERS.items()}

print(f"\nTree: {len(tips)} tips, {len(parent) - len(tips)} internal branches")
print("Kinetic origins:")
for name, u in origins.items():
    print(f"  {name:16s} branches: {', '.join(u['branches'])}")
print(f"Background: {len(bg_tips)} tip + {len(bg_internal)} internal branches")

# Matched null pseudo-origins: akinetic clades with the same tip count
bg_set = set(bg_tips)
cand_by_size = {}
for n in parent:
    cl = descendants_tips(n)
    if set(cl) <= bg_set:
        cand_by_size.setdefault(len(cl), []).append(tuple(sorted(cl)))
slots = []
for name, u in origins.items():
    c = cand_by_size.get(len(u["tips"]), [])
    if not c:
        print(f"  NOTE: no akinetic clade of {len(u['tips'])} tips to match "
              f"'{name}'; using single background tips instead.")
        c = cand_by_size[1]
    slots.append(c)
null_sets = set()
for combo in itertools.product(*slots):
    all_t = [t for cl in combo for t in cl]
    if len(all_t) == len(set(all_t)):
        null_sets.add(frozenset(combo))
null_sets = sorted(null_sets, key=lambda s: sorted(s))

# ─────────────────────────────────────────────────────────────────
# PER-RUN ANALYSIS
# ─────────────────────────────────────────────────────────────────
def unit_direction(u, change):
    """Per family: +1 if every tip's net change (stem -> tip) > 0,
    -1 if every one < 0, else 0. Also returns the net-change table."""
    net = pd.DataFrame({t.split("<")[0]: change[p].sum(axis=1)
                        for t, p in u["paths"].items()})
    d = pd.Series(0, index=change.index)
    d[(net > 0).all(axis=1)] = 1
    d[(net < 0).all(axis=1)] = -1
    return d, net


def unit_bp(u, change, bp):
    """Branch-probability flag per family: +1/-1 if the unit changes in a
    consistent direction AND one of its branches changing in that direction
    has branch probability < BP_THRESH; else 0. Also returns that minimum
    probability (NaN if CAFE5 did not score the branches)."""
    d, _ = unit_direction(u, change)
    cols = u["branches"]
    agree = change[cols].mul(d, axis=0) > 0
    pmin = bp[cols].where(agree).min(axis=1)
    flag = d.where(pmin < BP_THRESH, 0)
    return flag, pmin


def same_dir_matches(change, idx, cols, direction_series):
    ch = change.loc[idx, cols]
    same = ch.mul(direction_series, axis=0) > 0
    return same.apply(lambda r: [c.split("<")[0] if not c.startswith("<")
                                 else c for c in r.index[r]], axis=1)


def analyze(run):
    change, count, fam = run["change"], run["count"], run["fam"]
    res = pd.DataFrame(index=fam.index)
    res["family_pval"] = fam["pvalue"]
    res["family_sig"]  = fam["pvalue"] < P_THRESH

    for name, u in origins.items():
        d, net = unit_direction(u, change)
        res[f"{name}_dir"] = d.reindex(res.index)
        res[f"{name}_net"] = net.reindex(res.index).astype(int).astype(str) \
                                .agg("/".join, axis=1)
        for t in u["tips"]:
            if t in count.columns:
                res[f"{t.split('<')[0]}_count"] = count.loc[res.index, t]

    dir_cols = [f"{n}_dir" for n in origins]
    res["n_origins_expanded"]   = (res[dir_cols] > 0).sum(axis=1)
    res["n_origins_contracted"] = (res[dir_cols] < 0).sum(axis=1)
    res["shared_direction"] = 0
    res.loc[res["n_origins_expanded"] == K, "shared_direction"] = 1
    res.loc[res["n_origins_contracted"] == K, "shared_direction"] = -1

    sd = res["shared_direction"]
    res["bg_tip_matches"] = same_dir_matches(change, res.index, bg_tips, sd)
    res["bg_int_matches"] = same_dir_matches(change, res.index, bg_internal, sd)
    res["n_bg_tips_same"] = res["bg_tip_matches"].apply(len)
    res["n_bg_all_same"]  = res["n_bg_tips_same"] + res["bg_int_matches"].apply(len)
    res["frac_bg_tips_same"] = (res["n_bg_tips_same"] / len(bg_tips)).round(3)

    sister = pd.Series([[] for _ in res.index], index=res.index)
    for cols in sister_cols.values():
        m = same_dir_matches(change, res.index, cols, sd)
        sister = sister.combine(m, lambda a, b: a + b)
    res["sister_matches"] = sister.apply(lambda l: sorted(set(l)))
    res["unstable"] = res.index.isin(UNSTABLE_FAMILIES)

    shared = res["family_sig"] & (sd != 0)
    res["tier_all"]    = shared
    res["tier_strict"] = shared & (res["n_bg_tips_same"] == 0)
    res["tier_robust"] = shared & (res["n_bg_all_same"] == 0)
    res["tier_sister"] = shared & (res["sister_matches"].apply(len) == 0)
    res["tier_partial"] = res["family_sig"] & (
        ((res["n_origins_expanded"] == K - 1) & (res["n_origins_contracted"] == 0)) |
        ((res["n_origins_contracted"] == K - 1) & (res["n_origins_expanded"] == 0)))

    # Branch-probability screen (exploratory)
    bp = run["bp"]
    if bp is not None:
        for name, u in origins.items():
            f, pmin = unit_bp(u, change, bp)
            res[f"{name}_bp_min"] = pmin.reindex(res.index)
            res[f"{name}_bp_flag"] = f.reindex(res.index)
        bcols = [f"{n}_bp_flag" for n in origins]
        res["bp_shared_direction"] = 0
        res.loc[(res[bcols] > 0).all(axis=1), "bp_shared_direction"] = 1
        res.loc[(res[bcols] < 0).all(axis=1), "bp_shared_direction"] = -1
        bsd = res["bp_shared_direction"]
        # background tips with a low branch probability in the same direction
        bgp = bp.loc[res.index, bg_tips] < BP_THRESH
        bgs = change.loc[res.index, bg_tips].mul(bsd, axis=0) > 0
        both = bgp & bgs
        res["bp_bg_tip_matches"] = both.apply(
            lambda r: [c.split("<")[0] for c in r.index[r]], axis=1)
        res["n_bp_bg_tips_same"] = res["bp_bg_tip_matches"].apply(len)
        res["tier_bp_all"]    = res["family_sig"] & (bsd != 0)
        res["tier_bp_strict"] = res["tier_bp_all"] & (res["n_bp_bg_tips_same"] == 0)

    # Group-level test
    cache = {}

    def udir(tt):
        if tt not in cache:
            cache[tt] = unit_direction(make_unit(list(tt)), change)[0]
        return cache[tt]

    group = []
    null_rows = []
    for fam_label, fam_ids in [("significant", res.index[res["family_sig"]]),
                               ("all", res.index)]:
        obs = res.loc[fam_ids, dir_cols]
        obs_e = int((obs == 1).all(axis=1).sum())
        obs_c = int((obs == -1).all(axis=1).sum())
        null_e, null_c = [], []
        for s in null_sets:
            d = pd.concat([udir(cl) for cl in s], axis=1).loc[fam_ids]
            null_e.append(int((d == 1).all(axis=1).sum()))
            null_c.append(int((d == -1).all(axis=1).sum()))
        N = len(null_sets)
        for kind, nl in [("expansions", null_e), ("contractions", null_c)]:
            for k, v in enumerate(nl):
                null_rows.append({"run": run["dir"].name or str(run["dir"]),
                                  "score": "direction (change.tab)",
                                  "families": fam_label, "shared": kind,
                                  "null_set": k + 1, "count": v})
        for kind, o, nl in [("expansions", obs_e, null_e),
                            ("contractions", obs_c, null_c)]:
            s = pd.Series(nl)
            group.append({"families": f"{fam_label} (n={len(fam_ids)})",
                          "shared": kind, "observed": o,
                          "null_median": s.median(),
                          "null_range": f"{s.min()}-{s.max()}",
                          "p": (1 + (s >= o).sum()) / (1 + N)})
    group = pd.DataFrame(group)
    group.insert(0, "score", "direction (change.tab)")

    if bp is not None:
        bcache = {}

        def ubp(tt):
            if tt not in bcache:
                bcache[tt] = unit_bp(make_unit(list(tt)), change, bp)[0]
            return bcache[tt]

        bcols = [f"{n}_bp_flag" for n in origins]
        bgroup = []
        fam_ids = res.index[res["family_sig"]]
        obs = res.loc[fam_ids, bcols]
        obs_e = int((obs == 1).all(axis=1).sum())
        obs_c = int((obs == -1).all(axis=1).sum())
        null_e, null_c = [], []
        for s_ in null_sets:
            d = pd.concat([ubp(cl) for cl in s_], axis=1).loc[fam_ids]
            null_e.append(int((d == 1).all(axis=1).sum()))
            null_c.append(int((d == -1).all(axis=1).sum()))
        N = len(null_sets)
        for kind, nl in [("expansions", null_e), ("contractions", null_c)]:
            for k, v in enumerate(nl):
                null_rows.append({"run": run["dir"].name or str(run["dir"]),
                                  "score": f"branch prob < {BP_THRESH}",
                                  "families": "significant", "shared": kind,
                                  "null_set": k + 1, "count": v})
        for kind, o, nl in [("expansions", obs_e, null_e),
                            ("contractions", obs_c, null_c)]:
            s_ = pd.Series(nl)
            bgroup.append({"score": f"branch prob < {BP_THRESH}",
                           "families": f"significant (n={len(fam_ids)})",
                           "shared": kind, "observed": o,
                           "null_median": s_.median(),
                           "null_range": f"{s_.min()}-{s_.max()}",
                           "p": (1 + (s_ >= o).sum()) / (1 + N)})
        group = pd.concat([group, pd.DataFrame(bgroup)], ignore_index=True)
    group.attrs["null_rows"] = null_rows
    return res, group


results = []
for r in runs:
    res, group = analyze(r)
    results.append((r, res, group))

res, group = results[0][1], results[0][2]
TIERS = ["all", "strict", "robust", "sister", "partial"]
HAS_BP = "tier_bp_all" in res.columns
if HAS_BP:
    TIERS += ["bp_all", "bp_strict"]
n_runs = len(results)

# Replicate support: in how many runs is each family in each tier with the
# same shared direction as in the primary run?
if n_runs > 1:
    for t in TIERS:
        dcol = "bp_shared_direction" if t.startswith("bp_") else "shared_direction"
        hits = pd.Series(0, index=res.index)
        n_with = 0
        for _, rr, _ in results:
            if f"tier_{t}" not in rr or dcol not in rr:
                continue          # replicate without branch probabilities
            n_with += 1
            in_tier = rr[f"tier_{t}"].reindex(res.index, fill_value=False)
            same_dir = rr[dcol].reindex(res.index) == res[dcol]
            hits += (in_tier & same_dir).astype(int)
        res[f"reps_{t}"] = hits.astype(str) + f"/{n_with}"

# ─────────────────────────────────────────────────────────────────
# OUTPUT
# ─────────────────────────────────────────────────────────────────
def fmt(df):
    out = df.copy()
    for c in ["bg_tip_matches", "bg_int_matches", "sister_matches", "bp_bg_tip_matches"]:
        if c in out:
            out[c] = out[c].apply(lambda l: ",".join(l) if l else "-")
    out.insert(0, "direction", out["shared_direction"]
               .map({1: "expansion", -1: "contraction", 0: "mixed"}))
    return out.sort_values("family_pval")


SHOW = ["direction", "family_pval"] + [f"{n}_net" for n in origins] + \
       ["n_bg_tips_same", "sister_matches", "unstable"]

print("\n" + "=" * 70)
print(f"  PER-FAMILY RESULTS (primary run): family-wide p < {P_THRESH} AND "
      f"all {K} kinetic origins change in the same direction")
if n_runs > 1:
    print(f"  'reps' = number of the {n_runs} replicates recovering the "
          "family in that tier, same direction")
print("=" * 70)
for label, key, desc in [
        ("ALL SHARED", "all", "no background filter"),
        ("STRICT", "strict", "no background tip changes in the same direction"),
        ("ROBUST", "robust", "no background tip OR internal branch in the "
                             "same direction"),
        ("SISTER-PAIR", "sister", "no akinetic sister changes in the "
                                  "same direction")]:
    df = res[res[f"tier_{key}"]]
    n_e = (df["shared_direction"] == 1).sum()
    n_c = (df["shared_direction"] == -1).sum()
    extra = ""
    if n_runs > 1 and len(df):
        n_all = (df[f"reps_{key}"] == f"{n_runs}/{n_runs}").sum()
        extra = f"  — {n_all} recovered in all {n_runs} replicates"
    print(f"\n── {label} ({desc}): {len(df)}  "
          f"[{n_e} expansions, {n_c} contractions]{extra}")
    if len(df):
        cols = SHOW + ([f"reps_{key}"] if n_runs > 1 else [])
        print(fmt(df)[cols].to_string())
        fmt(df).to_csv(OUT_DIR / f"kinetic_shared_{key}.tsv", sep="\t",
                       index_label="FamilyID")

n_partial = res["tier_partial"].sum()
print(f"\n── PARTIAL ({K - 1} of {K} origins same direction, none opposite): "
      f"{n_partial} families (tier_partial in kinetic_all_families.tsv)")

fmt(res).to_csv(OUT_DIR / "kinetic_all_families.tsv", sep="\t",
                index_label="FamilyID")

if HAS_BP:
    print("\n" + "=" * 70)
    print(f"  BRANCH-PROBABILITY SCREEN (exploratory; branch prob < {BP_THRESH})")
    print("  Branch probabilities are used as a score, NOT as p-values. Families")
    print("  listed are candidates for follow-up, not lineage-significant results.")
    print("=" * 70)
    BP_SHOW = ["direction_bp", "family_pval"] + \
              [f"{n}_bp_min" for n in origins] + [f"{n}_net" for n in origins] + \
              ["n_bp_bg_tips_same", "tier_all"]
    for label, key, desc in [
            ("BP ALL", "bp_all", "every origin flagged, same direction"),
            ("BP STRICT", "bp_strict", "... and no background tip flagged "
                                       "in that direction")]:
        df = res[res[f"tier_{key}"]]
        n_e = (df["bp_shared_direction"] == 1).sum()
        n_c = (df["bp_shared_direction"] == -1).sum()
        extra = ""
        if n_runs > 1 and len(df):
            n_bp = int(df[f"reps_{key}"].iloc[0].split("/")[1])
            n_rep = (df[f"reps_{key}"] == f"{n_bp}/{n_bp}").sum()
            extra = (f"  — {n_rep} recovered in all {n_bp} replicates with "
                     "branch probabilities")
        print(f"\n── {label} ({desc}): {len(df)}  "
              f"[{n_e} expansions, {n_c} contractions]{extra}")
        if len(df):
            out = fmt(df)
            out.insert(0, "direction_bp", out["bp_shared_direction"]
                       .map({1: "expansion", -1: "contraction", 0: "mixed"}))
            cols = BP_SHOW + ([f"reps_{key}"] if n_runs > 1 else [])
            print(out[cols].to_string())
            out.to_csv(OUT_DIR / f"kinetic_shared_{key}.tsv", sep="\t",
                       index_label="FamilyID")
    n_scored = res.loc[res["family_sig"], [f"{n}_bp_min" for n in origins]] \
        .notna().any(axis=1).sum()
    print(f"\n  Significant families with branch probabilities on kinetic "
          f"branches: {n_scored} / {int(res['family_sig'].sum())}")

shared = res[res["tier_all"]]
tally = Counter(t for l in shared["bg_tip_matches"] for t in l)
if tally:
    print("\nBackground tips most often sharing the kinetic direction "
          f"(out of {len(shared)} shared families):")
    for t, n in tally.most_common():
        print(f"  {t:32s} {n}")

print("\n" + "=" * 70)
print("  GROUP-LEVEL TEST: do kinetic origins share expansions/contractions")
print("  more often than structurally matched sets of akinetic lineages?")
print("=" * 70)
print(f"  {len(null_sets)} matched null sets; smallest attainable p = "
      f"{1 / (1 + len(null_sets)):.4f}\n")
all_groups = []
for r, _, g in results:
    g = g.copy()
    g.insert(0, "run", r["dir"].name or str(r["dir"]))
    all_groups.append(g)
all_groups = pd.concat(all_groups, ignore_index=True)
if n_runs > 1:
    print("  (first run listed = primary)")
print(all_groups.to_string(index=False))
all_groups.to_csv(OUT_DIR / "kinetic_group_test.tsv", sep="\t", index=False)
# per-null-set counts (for plotting the null distributions)
pd.DataFrame([row for _, _, g in results for row in g.attrs.get("null_rows", [])]) \
    .to_csv(OUT_DIR / "kinetic_null_counts.tsv", sep="\t", index=False)

pd.DataFrame({"null_set": [" | ".join("+".join(t.split("<")[0] for t in cl)
                                      for cl in sorted(s)) for s in null_sets]}) \
    .to_csv(OUT_DIR / "kinetic_null_sets.tsv", sep="\t", index=False)
print("\nOutputs: kinetic_all_families.tsv, kinetic_shared_<tier>.tsv, "
      "kinetic_group_test.tsv, kinetic_null_counts.tsv, kinetic_null_sets.tsv"
      + (", kinetic_replicates.tsv" if n_runs > 1 else ""))