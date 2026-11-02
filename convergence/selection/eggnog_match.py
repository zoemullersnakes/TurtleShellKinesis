# Add eggNOG annotations to the HyPhy summary
#
# The HyPhy 'reference_header' (e.g. Chrysemys_picta_bellii-LOC101949731) is the
# same ID as eggNOG's 'query', so the two tables are joined on that column.
# Every column of the HyPhy sheet is kept, including any names you entered by
# hand in 'gene_name'; the eggNOG columns are added after 'reference_header'.
#
# Added columns
#   ncbi_symbol    symbol from the Chrysemys RefSeq header (blank for LOC genes)
#   eggnog_name    eggNOG Preferred_name (blank where eggNOG gave '-')
#   names_agree    True/False where both exist, so you can review disagreements
#   final_symbol   NCBI symbol if there is one, otherwise the eggNOG name,
#                  otherwise the LOC ID
#   symbol_source  which of the three final_symbol came from
#   + eggNOG COG_category, GOs, KEGG_ko, KEGG_Pathway, PFAMs

import pandas as pd


def eggnog_match(eggnog_xlsx, hyphy_xlsx, out, sheet="Summary"):
    # Step 1: read eggNOG, keeping the useful columns
    eggnog = pd.read_excel(eggnog_xlsx)
    keep = ["query", "Preferred_name", "COG_category", "GOs",
            "KEGG_ko", "KEGG_Pathway", "PFAMs"]
    egg = eggnog[[c for c in keep if c in eggnog.columns]].copy()
    egg = egg.rename(columns={"Preferred_name": "eggnog_name"})
    egg["eggnog_name"] = egg["eggnog_name"].replace("-", pd.NA)

    # Step 2: read the HyPhy summary
    hyphy = pd.read_excel(hyphy_xlsx, sheet_name=sheet)

    # Step 3: join on the Chrysemys protein ID
    merged = hyphy.merge(egg, left_on="reference_header", right_on="query",
                         how="left").drop(columns="query")

    # Step 4: NCBI symbol vs eggNOG name, and one final symbol per gene
    sym = merged["reference_header"].str.split("-", n=1).str[1]
    is_loc = sym.str.match(r"LOC\d+$", na=False)
    merged["ncbi_symbol"] = sym.where(~is_loc)
    both = merged["ncbi_symbol"].notna() & merged["eggnog_name"].notna()
    merged["names_agree"] = pd.NA
    merged.loc[both, "names_agree"] = (
        merged.loc[both, "ncbi_symbol"].str.upper()
        == merged.loc[both, "eggnog_name"].str.upper())

    egg_ok = merged["eggnog_name"].notna() & ~merged["eggnog_name"].astype(str).str.match(r"LOC\d+$")
    merged["final_symbol"] = merged["ncbi_symbol"]
    merged["symbol_source"] = "NCBI"
    # names you already entered by hand in gene_name take priority over eggNOG
    if "gene_name" in merged.columns:
        manual = merged["gene_name"].astype(str)
        use_manual = (merged["ncbi_symbol"].isna() & merged["gene_name"].notna()
                      & ~manual.str.match(r"LOC\d+$"))
        merged.loc[use_manual, "final_symbol"] = merged.loc[use_manual, "gene_name"]
        merged.loc[use_manual, "symbol_source"] = "manual"
    use_egg = merged["final_symbol"].isna() & egg_ok
    merged.loc[use_egg, "final_symbol"] = merged.loc[use_egg, "eggnog_name"]
    merged.loc[use_egg, "symbol_source"] = "eggNOG"
    still = merged["final_symbol"].isna()
    merged.loc[still, "final_symbol"] = sym[still]
    merged.loc[still, "symbol_source"] = "LOC (unnamed)"

    # Step 5: put the new columns right after reference_header
    new = ["ncbi_symbol", "eggnog_name", "names_agree", "final_symbol", "symbol_source"]
    rest = [c for c in merged.columns if c not in new]
    at = rest.index("reference_header") + 1
    merged = merged[rest[:at] + new + rest[at:]]

    # Step 6: write out
    if out.endswith(".xlsx"):
        merged.to_excel(out, index=False, sheet_name=sheet)
    else:
        merged.to_csv(out, sep="\t", index=False)

    # Summary of what happened
    n = len(merged)
    print(f"HyPhy rows: {n}; matched to eggNOG: {merged['reference_header'].isin(eggnog['query']).sum()}")
    print(merged["symbol_source"].value_counts().to_string())
    print(f"NCBI and eggNOG names agree for {int(merged['names_agree'].eq(True).sum())} "
          f"of {int(both.sum())} genes that have both")
    return merged


eggnog_xlsx = "/Users/zlm9173/Desktop/Turtles/emapper_annotations.xlsx"
hyphy_xlsx = "/Users/zlm9173/Desktop/Turtles/hyphy_out.xlsx"
out = "/Users/zlm9173/Desktop/Turtles/hyphy_out_eggnog.xlsx"

if __name__ == "__main__":
    eggnog_match(eggnog_xlsx, hyphy_xlsx, out)