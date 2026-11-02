#!/usr/bin/env python3
"""
hyphy_to_xlsx.py
================
Collects HyPhy RELAX, aBSREL and BUSTED JSON outputs into one Excel workbook.

Sheets
  Summary          one row per ortholog: RELAX k / p / q and direction,
                   BUSTED p / q, aBSREL significant branches; sig_* flags
  RELAX            per-gene test results (+ BH q-value across genes)
  RELAX_branches   per-branch k (general descriptive model) with Test /
                   Reference label
  BUSTED           per-gene test results, omega classes of the Test branches
                   (+ BH q-value)
  BUSTED-PH        per-gene: test / background / shared-distribution LRTs,
                   BH q for test and difference, omega classes of each set,
                   interpretation category, ph_trait_associated flag
                   (test q < alpha, background raw p >= alpha, difference
                   q < alpha, larger max omega on test, test omega < --omega-cap)
  aBSREL           per-gene summary (branches tested / significant)
  aBSREL_branches  per-branch LRT, uncorrected and Holm-corrected p (HyPhy's
                   own correction across branches), omega classes
  Files            every JSON read, the ortholog name assigned, and any
                   problem (unreadable, missing sections)

Every gene-level sheet reports how many branches were in the Test (foreground)
set, so a run with the wrong foreground (e.g. every branch tested) is visible.

USAGE
  python hyphy_to_xlsx.py --relax relax_dir/ --absrel absrel_dir/ \
                          --busted busted_dir/ --bustedph bustedph_dir/ \
                          -o hyphy_results.xlsx
  python hyphy_to_xlsx.py --all results_dir/ -o hyphy_results.xlsx
      (--all scans one folder for *.relax.json, *.absrel.json, *.busted.json,
       *.bustedph.json)

Ortholog name = file name up to the first "_" or "." (OG0001234_cds.RELAX.json
-> OG0001234); change with --name-regex.

Requires pandas and openpyxl (pip install --user pandas openpyxl).
"""

import argparse
import json
import os
import re
import sys

import pandas as pd

SUFFIX = {"RELAX": ".relax.json", "aBSREL": ".absrel.json", "BUSTED": ".busted.json", "BUSTED-PH": ".bustedph.json"}


# ─────────────────────────── helpers ───────────────────────────
def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values; NaN entries stay NaN."""
    s = pd.Series(pvals, dtype="float64")
    ok = s.dropna().sort_values()
    m = len(ok)
    if m == 0:
        return s
    ranked = ok * m / pd.Series(range(1, m + 1), index=ok.index)
    q = ranked[::-1].cummin()[::-1].clip(upper=1.0)
    return q.reindex(s.index)


def find_files(folder, method):
    if not folder:
        return []
    if not os.path.isdir(folder):
        sys.exit(f"ERROR: not a folder: {folder}")
    suf = SUFFIX[method]
    return sorted(os.path.join(folder, f) for f in os.listdir(folder)
                  if f.lower().endswith(suf))


def load(path, log, method, name_re):
    fname = os.path.basename(path)
    m = re.match(name_re, fname)
    ortholog = m.group(1) if m else fname
    entry = {"method": method, "file": fname, "ortholog": ortholog, "problem": ""}
    log.append(entry)
    try:
        with open(path) as fh:
            data = json.load(fh)
    except Exception as e:  # truncated / still-running jobs
        entry["problem"] = f"unreadable JSON: {e}"
        return ortholog, None, entry
    return ortholog, data, entry


def tested_labels(data):
    """{branch: label} from the 'tested' section (labels lower-cased)."""
    t = data.get("tested", {})
    t = t.get("0", t) if isinstance(t, dict) else {}
    return {b: str(v).lower() for b, v in t.items()} if isinstance(t, dict) else {}


def omega_classes(dist):
    """Rate-distribution dict {"0": {"omega", "proportion"}, ...} or list of
    [omega, weight] pairs -> (text, max omega, weight of max class)."""
    pairs = []
    if isinstance(dist, dict):
        for k in sorted(dist, key=lambda x: int(x) if str(x).isdigit() else 0):
            v = dist[k]
            if isinstance(v, dict) and "omega" in v:
                pairs.append((v.get("omega"), v.get("proportion")))
    elif isinstance(dist, list):
        for v in dist:
            if isinstance(v, (list, tuple)) and len(v) >= 2:
                pairs.append((v[0], v[1]))
    pairs = [(w, p) for w, p in pairs if w is not None]
    if not pairs:
        return "", None, None
    text = "; ".join(f"ω={w:.3g} ({100 * p:.1f}%)" if p is not None else f"ω={w:.3g}"
                     for w, p in pairs)
    w_max, p_max = max(pairs, key=lambda x: x[0])
    return text, w_max, p_max


def get_fit_dist(data, fit_names, part_names):
    fits = data.get("fits", {})
    for fn in fit_names:
        rd = fits.get(fn, {}).get("Rate Distributions", {})
        for pn in part_names:
            if pn in rd:
                return rd[pn]
    return None


# ─────────────────────────── gene names ───────────────────────────
DEFAULT_GENE_PATTERNS = [
    r"\[gene=([^\]]+)\]",            # NCBI style: [gene=TP53]
    r"\bgene[=:]\s*([^\s;|\]]+)",    # gene=TP53 / gene:TP53
    r"\bGN=([^\s;|]+)",              # UniProt style
]
FASTA_EXT = (".fa", ".fasta", ".faa", ".fna", ".fas", ".ffn", ".aln")


def load_gene_names(folder, species, gene_regex, name_re):
    """{ortholog: (gene_name, header)} from the header of the `species`
    sequence in each original ortholog FASTA."""
    if not os.path.isdir(folder):
        sys.exit(f"ERROR: not a folder: {folder}")
    patterns = [gene_regex] if gene_regex else DEFAULT_GENE_PATTERNS
    sp = species.lower()
    out = {}
    for f in sorted(os.listdir(folder)):
        if not f.lower().endswith(FASTA_EXT):
            continue
        m = re.match(name_re, f)
        og = m.group(1) if m else f
        header = None
        with open(os.path.join(folder, f)) as fh:
            for line in fh:
                if line.startswith(">") and sp in line.lower():
                    header = line[1:].strip()
                    break
        if header is None:
            out[og] = ("", "")
            continue
        gene = ""
        for pat in patterns:
            g = re.search(pat, header)
            if g:
                gene = g.group(1)
                break
        out[og] = (gene, header)
    return out


def add_gene_column(df, genes, with_header=False):
    if df is None or df.empty or "ortholog" not in df:
        return df
    pos = df.columns.get_loc("ortholog") + 1
    df.insert(pos, "gene_name", df["ortholog"].map(lambda o: genes.get(o, ("", ""))[0]))
    if with_header:
        df.insert(pos + 1, "reference_header",
                  df["ortholog"].map(lambda o: genes.get(o, ("", ""))[1]))
    return df


# ─────────────────────────── parsers ───────────────────────────
def parse_relax(files, log, name_re):
    genes, branches = [], []
    for path in files:
        og, data, entry = load(path, log, "RELAX", name_re)
        if data is None:
            continue
        tr = data.get("test results", {})
        labels = tested_labels(data)
        k = tr.get("relaxation or intensification parameter")
        genes.append({
            "ortholog": og,
            "relax_k": k,
            "relax_LRT": tr.get("LRT"),
            "relax_p": tr.get("p-value"),
            "relax_direction": ("intensified" if k is not None and k > 1 else
                                "relaxed" if k is not None and k < 1 else ""),
            "relax_n_test_branches": sum(v == "test" for v in labels.values()),
            "relax_n_reference_branches": sum(v == "reference" for v in labels.values()),
            "relax_test_omegas": omega_classes(get_fit_dist(
                data, ["RELAX alternative", "Alternative"], ["Test"]))[0],
            "relax_reference_omegas": omega_classes(get_fit_dist(
                data, ["RELAX alternative", "Alternative"], ["Reference"]))[0],
            "file": entry["file"],
        })
        if not tr:
            entry["problem"] = "no 'test results' section"
        ba = data.get("branch attributes", {}).get("0", {})
        n_k = 0
        for b, attrs in ba.items():
            if not isinstance(attrs, dict):
                continue
            bk = attrs.get("k (general descriptive)")
            if bk is None:
                continue
            n_k += 1
            branches.append({"ortholog": og, "branch": b,
                             "branch_set": labels.get(b, ""),
                             "branch_k": bk, "gene_k": k,
                             "gene_p": tr.get("p-value")})
        if n_k == 0 and not entry["problem"]:
            entry["problem"] = "no per-branch k (general descriptive model not run?)"
    g = pd.DataFrame(genes)
    if len(g):
        g.insert(g.columns.get_loc("relax_p") + 1, "relax_q", bh_qvalues(g["relax_p"]))
    return g, pd.DataFrame(branches)


def parse_busted(files, log, name_re):
    genes = []
    for path in files:
        og, data, entry = load(path, log, "BUSTED", name_re)
        if data is None:
            continue
        tr = data.get("test results", {})
        labels = tested_labels(data)
        txt, wmax, pmax = omega_classes(get_fit_dist(
            data, ["Unconstrained model", "Unconstrained"], ["Test"]))
        bg_txt = omega_classes(get_fit_dist(
            data, ["Unconstrained model", "Unconstrained"], ["Background"]))[0]
        genes.append({
            "ortholog": og,
            "busted_LRT": tr.get("LRT"),
            "busted_p": tr.get("p-value"),
            "busted_test_max_omega": wmax,
            "busted_test_max_omega_prop": pmax,
            "busted_test_omegas": txt,
            "busted_background_omegas": bg_txt,
            "busted_n_test_branches": sum(v == "test" for v in labels.values()),
            "busted_n_background_branches": sum(v == "background" for v in labels.values()),
            "file": entry["file"],
        })
        if not tr:
            entry["problem"] = "no 'test results' section"
    g = pd.DataFrame(genes)
    if len(g):
        g.insert(g.columns.get_loc("busted_p") + 1, "busted_q", bh_qvalues(g["busted_p"]))
    return g


def parse_absrel(files, log, name_re, alpha):
    genes, branches = [], []
    for path in files:
        og, data, entry = load(path, log, "aBSREL", name_re)
        if data is None:
            continue
        tr = data.get("test results", {})
        labels = tested_labels(data)
        ba = data.get("branch attributes", {}).get("0", {})
        sig = []
        n_rows = 0
        for b, a in ba.items():
            if not isinstance(a, dict) or "Corrected P-value" not in a:
                continue          # untested branches carry no test results
            n_rows += 1
            txt, wmax, pmax = omega_classes(a.get("Rate Distributions"))
            pc = a.get("Corrected P-value")
            branches.append({
                "ortholog": og, "branch": b,
                "branch_set": labels.get(b, ""),
                "LRT": a.get("LRT"),
                "p_uncorrected": a.get("Uncorrected P-value"),
                "p_holm": pc,
                "significant": pc is not None and pc < alpha,
                "rate_classes": a.get("Rate classes"),
                "max_omega": wmax, "max_omega_prop": pmax,
                "omega_classes": txt,
                "baseline_omega": a.get("Baseline MG94xREV omega ratio"),
            })
            if pc is not None and pc < alpha:
                sig.append(b)
        genes.append({
            "ortholog": og,
            "absrel_n_tested": tr.get("tested", n_rows),
            "absrel_n_test_label": sum(v == "test" for v in labels.values()),
            "absrel_n_significant": len(sig),
            "absrel_significant_branches": ",".join(sig),
            "absrel_min_p_holm": min((r["p_holm"] for r in branches
                                      if r["ortholog"] == og and r["p_holm"] is not None),
                                     default=None),
            "file": entry["file"],
        })
        if n_rows == 0:
            entry["problem"] = "no tested branches with results"
    return pd.DataFrame(genes), pd.DataFrame(branches)


def _find_test(data, kind):
    """LRT / p-value dict for one BUSTED-PH test, from either output format.
    kind: 'fg' | 'bg' | 'diff'.
      built-in `hyphy busted-ph`: 'Background selection test results',
        'Comparative selection test results', foreground as 'test results'
        or 'Foreground selection test results'
      standalone BUSTED-PH.bf:    'test results', 'test results background',
        'test results shared distributions'"""
    words = {"bg": ("background",), "diff": ("comparative", "shared")}
    for k, v in data.items():
        kl = k.lower()
        if "test results" not in kl or not isinstance(v, dict):
            continue
        is_bg = "background" in kl
        is_diff = any(w in kl for w in words["diff"])
        if kind == "bg" and is_bg and not is_diff:
            return v
        if kind == "diff" and is_diff:
            return v
        if kind == "fg" and not is_bg and not is_diff:
            return v
    return None


def _ph_block_p(data, kind):
    """p-value from the built-in version's 'BUSTED-PH' summary block."""
    blk = data.get("BUSTED-PH", {})
    up = blk.get("Uncorrected P-value", {}) if isinstance(blk, dict) else {}
    key = {"fg": "FG", "bg": "BG", "diff": "Comparative"}[kind]
    return up.get(key) if isinstance(up, dict) else None


def _ph_category(sig_t, sig_b, sig_d):
    if not sig_t:
        return "no selection on test"
    if sig_b and sig_d:
        return "selection in both, different regimes"
    if sig_b:
        return "selection in both, similar regimes"
    if sig_d:
        return "trait-associated"
    return "test only, regimes not different"


def parse_bustedph(files, log, name_re, alpha, omega_cap):
    genes = []
    for path in files:
        og, data, entry = load(path, log, "BUSTED-PH", name_re)
        if data is None:
            continue
        res = {}
        for kind in ("fg", "bg", "diff"):
            tr = _find_test(data, kind) or {}
            p = tr.get("p-value")
            if p is None:
                p = _ph_block_p(data, kind)
            res[kind] = (tr.get("LRT"), p)
        blk = data.get("BUSTED-PH", {}) if isinstance(data.get("BUSTED-PH"), dict) else {}
        labels = tested_labels(data)
        fits = ["Unconstrained model", "Unconstrained"]
        t_txt, t_w, t_p = omega_classes(get_fit_dist(data, fits, ["Test", "Foreground", "FG"]))
        b_txt, b_w, b_p = omega_classes(get_fit_dist(data, fits, ["Background", "BG"]))
        genes.append({
            "ortholog": og,
            "ph_test_LRT": res["fg"][0],
            "ph_test_p": res["fg"][1],
            "ph_bg_LRT": res["bg"][0],
            "ph_bg_p": res["bg"][1],
            # background test is skipped when the unconstrained fit has no
            # omega>1 class on background: no evidence of background selection
            "ph_bg_tested": res["bg"][1] is not None,
            "ph_diff_LRT": res["diff"][0],
            "ph_diff_p": res["diff"][1],
            "ph_hyphy_verdict": ("associated" if "**no**" not in str(blk.get("Summary", ""))
                                 and blk.get("Summary") else
                                 "not associated" if blk.get("Summary") else ""),
            "ph_hyphy_summary": blk.get("Summary", ""),
            "ph_test_max_omega": t_w, "ph_test_max_omega_prop": t_p,
            "ph_bg_max_omega": b_w, "ph_bg_max_omega_prop": b_p,
            "ph_test_omegas": t_txt, "ph_background_omegas": b_txt,
            "ph_n_test_branches": sum(v == "test" for v in labels.values()),
            "ph_n_background_branches": sum(v == "background" for v in labels.values()),
            "ph_n_nuisance_branches": sum(v == "nuisance" for v in labels.values()),
            "file": entry["file"],
        })
        if res["fg"][1] is None:
            entry["problem"] = "no foreground test result"
        elif res["diff"][1] is None:
            entry["problem"] = "no comparative / shared-distribution test result"
    g = pd.DataFrame(genes)
    if not len(g):
        return g
    for col in ("test", "diff"):
        g.insert(g.columns.get_loc(f"ph_{col}_p") + 1, f"ph_{col}_q",
                 bh_qvalues(g[f"ph_{col}_p"]))
    # Test-branch selection and the difference are the discoveries, so they use
    # BH q across genes. Background selection is what rules a gene OUT, so it
    # uses the raw p: stricter for calling a gene trait-associated.
    sig_t = g["ph_test_q"] < alpha
    sig_b = g["ph_bg_p"].fillna(1) < alpha
    sig_d = g["ph_diff_q"] < alpha
    g["ph_category"] = [_ph_category(*x) for x in zip(sig_t, sig_b, sig_d)]
    tw = pd.to_numeric(g["ph_test_max_omega"], errors="coerce")
    bw = pd.to_numeric(g["ph_bg_max_omega"], errors="coerce").fillna(1)
    # NA when the rate classes could not be read, so they don't block the flag
    g["ph_test_stronger"] = (tw > bw).astype("boolean").mask(tw.isna())
    g["ph_omega_outlier"] = tw >= omega_cap
    g["ph_trait_associated"] = ((g["ph_category"] == "trait-associated")
                                & g["ph_test_stronger"].fillna(True).astype(bool)
                                & ~g["ph_omega_outlier"])
    if tw.isna().all():
        print("   ⚠️  BUSTED-PH: no omega rate classes found under fits/Unconstrained model/"
              "Rate Distributions; direction and omega-outlier checks were skipped")
    return g


# ─────────────────────────── summary ───────────────────────────
def build_summary(relax, busted, absrel, alpha, bustedph=None):
    parts = []
    if bustedph is not None and len(bustedph):
        parts.append(bustedph[["ortholog", "ph_category", "ph_trait_associated", "ph_hyphy_verdict",
                               "ph_test_q", "ph_bg_p", "ph_diff_q",
                               "ph_test_max_omega", "ph_bg_max_omega", "ph_omega_outlier"]])
    if len(relax):
        parts.append(relax[["ortholog", "relax_k", "relax_direction", "relax_p", "relax_q"]])
    if len(busted):
        parts.append(busted[["ortholog", "busted_p", "busted_q", "busted_test_max_omega",
                             "busted_test_max_omega_prop"]])
    if len(absrel):
        parts.append(absrel[["ortholog", "absrel_n_significant",
                             "absrel_significant_branches", "absrel_min_p_holm"]])
    if not parts:
        return pd.DataFrame()
    for p in parts:   # one row per ortholog per method
        dup = p["ortholog"].duplicated(keep=False)
        if dup.any():
            print(f"⚠️  duplicate ortholog names within one method: "
                  f"{sorted(p.loc[dup, 'ortholog'].unique())[:5]} (check --name-regex)")
    s = parts[0].drop_duplicates("ortholog")
    for p in parts[1:]:
        s = s.merge(p.drop_duplicates("ortholog"), on="ortholog", how="outer")
    if "relax_q" in s:
        s["sig_RELAX"] = s["relax_q"] < alpha
    if "busted_q" in s:
        s["sig_BUSTED"] = s["busted_q"] < alpha
    if "absrel_n_significant" in s:
        s["sig_aBSREL"] = s["absrel_n_significant"].fillna(0) > 0
    if "ph_trait_associated" in s:
        s["sig_BUSTEDPH"] = s["ph_trait_associated"].fillna(False).astype(bool)
    flags = [c for c in s.columns if c.startswith("sig_")]
    s["n_methods_sig"] = s[flags].fillna(False).astype(bool).sum(axis=1)
    sort_cols = ["n_methods_sig"] + [c for c in ["ph_test_q", "busted_q", "relax_q"] if c in s]
    return s.sort_values(sort_cols, ascending=[False] + [True] * (len(sort_cols) - 1))


# ─────────────────────────── writer ───────────────────────────
def write_xlsx(sheets, out):
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        for name, df in sheets.items():
            if df is None or df.empty:
                continue
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            ws.freeze_panes = "B2"
            ws.auto_filter.ref = ws.dimensions
            hdr_fill = PatternFill("solid", fgColor="DDE4EE")
            for c in ws[1]:
                c.font = Font(bold=True)
                c.fill = hdr_fill
            for i, col in enumerate(df.columns, start=1):
                letter = get_column_letter(i)
                width = max(len(str(col)), *(len(str(v)) for v in df[col].head(500))) + 2
                ws.column_dimensions[letter].width = min(width, 60)
                low = col.lower()
                fmt = None
                if low.endswith("_p") or low.endswith("_q") or "p_" in low or low.startswith("p_") \
                        or low in ("gene_p",) or "p_holm" in low:
                    fmt = "0.000E+00"
                elif any(t in low for t in ("lrt", "omega", "_k", "prop")) and "omegas" not in low \
                        and "classes" not in low:
                    fmt = "0.0000"
                if fmt:
                    for cell in ws.iter_rows(min_row=2, min_col=i, max_col=i):
                        cell[0].number_format = fmt


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--relax", help="folder of *.RELAX.json")
    ap.add_argument("--absrel", help="folder of *.ABSREL.json")
    ap.add_argument("--busted", help="folder of *.BUSTED.json")
    ap.add_argument("--bustedph", help="folder of *.bustedph.json")
    ap.add_argument("--omega-cap", type=float, default=1e4,
                    help="BUSTED-PH genes with test omega >= this are flagged as likely "
                         "alignment error and excluded from ph_trait_associated (default 1e4)")
    ap.add_argument("--all", help="one folder containing all kinds")
    ap.add_argument("-o", "--output", required=True, help="output .xlsx")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="significance threshold for q-values and aBSREL corrected p (default 0.05)")
    ap.add_argument("--name-regex", default=r"^([^_.]+)",
                    help="regex with one group giving the ortholog name from the file name")
    ap.add_argument("--genes", help="folder of the original ortholog FASTA files; "
                    "adds a gene_name column taken from the reference species' header")
    ap.add_argument("--gene-species", default="Chrysemys",
                    help="text identifying the reference species' header (default Chrysemys)")
    ap.add_argument("--gene-regex",
                    help="regex with one group extracting the gene name from that header "
                         "(default tries [gene=X], gene=X, GN=X)")
    a = ap.parse_args()

    if a.all:
        a.relax = a.relax or a.all
        a.absrel = a.absrel or a.all
        a.busted = a.busted or a.all
        a.bustedph = a.bustedph or a.all
    if not (a.relax or a.absrel or a.busted or a.bustedph):
        ap.error("give at least one of --relax / --absrel / --busted / --bustedph / --all")
    if not a.output.lower().endswith(".xlsx"):
        a.output += ".xlsx"

    log = []
    files = {m: find_files(d, m) for m, d in
             [("RELAX", a.relax), ("aBSREL", a.absrel), ("BUSTED", a.busted),
              ("BUSTED-PH", a.bustedph)]}
    for m, f in files.items():
        print(f"🔍 {m}: {len(f)} JSON files")

    relax, relax_br = parse_relax(files["RELAX"], log, a.name_regex)
    busted = parse_busted(files["BUSTED"], log, a.name_regex)
    absrel, absrel_br = parse_absrel(files["aBSREL"], log, a.name_regex, a.alpha)
    bustedph = parse_bustedph(files["BUSTED-PH"], log, a.name_regex, a.alpha, a.omega_cap)
    summary = build_summary(relax, busted, absrel, a.alpha, bustedph)
    if len(bustedph):
        print("   BUSTED-PH categories: " + ", ".join(
            f"{k} ({v})" for k, v in bustedph["ph_category"].value_counts().items()))
        print(f"   BUSTED-PH background test skipped (no ω>1 on background): "
              f"{int((~bustedph['ph_bg_tested']).sum())} genes")

    # Foreground sanity check: the number of Test branches should be the same
    # for every gene with the full taxon set; a jump suggests the wrong flags.
    for df, col, m in [(relax, "relax_n_test_branches", "RELAX"),
                       (busted, "busted_n_test_branches", "BUSTED"),
                       (bustedph, "ph_n_test_branches", "BUSTED-PH"),
                       (absrel, "absrel_n_tested", "aBSREL")]:
        if len(df) and col in df:
            vc = df[col].value_counts().sort_index()
            print(f"   {m} branches tested per gene: " +
                  ", ".join(f"{k} ({v} genes)" for k, v in vc.items()))

    if a.genes:
        genes = load_gene_names(a.genes, a.gene_species, a.gene_regex, a.name_regex)
        n_hdr = sum(1 for g, h in genes.values() if h)
        n_gene = sum(1 for g, h in genes.values() if g)
        print(f"🧬 Gene names: {len(genes)} ortholog files, {n_hdr} with a "
              f"{a.gene_species} sequence, {n_gene} with a gene name parsed")
        if n_hdr and not n_gene:
            ex = next(h for g, h in genes.values() if h)
            print(f"   ⚠️  no gene names parsed; example header: {ex}\n"
                  f"      set --gene-regex to match it")
        summary = add_gene_column(summary, genes, with_header=True)
        for df in (relax, relax_br, busted, bustedph, absrel, absrel_br):
            add_gene_column(df, genes)

    write_xlsx({"Summary": summary, "RELAX": relax, "RELAX_branches": relax_br,
                "BUSTED": busted, "BUSTED-PH": bustedph, "aBSREL": absrel, "aBSREL_branches": absrel_br,
                "Files": pd.DataFrame(log)}, a.output)

    bad = [e for e in log if e["problem"]]
    print(f"✅ Wrote {a.output}: {len(summary)} orthologs"
          + (f"; ⚠️ {len(bad)} file(s) with problems (see 'Files' sheet)" if bad else ""))
    for flag in [c for c in summary.columns if c.startswith("sig_")]:
        print(f"   {flag}: {int(summary[flag].fillna(False).astype(bool).sum())}")


if __name__ == "__main__":
    main()