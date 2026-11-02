import pandas as pd
import os
import shutil

# Configuration
ORTHOGROUPS_TSV = "Orthogroups/Orthogroups.tsv"
#ORTHOGROUP_SEQS_DIR = "OrthoFinder/Results/Orthogroup_Sequences"  # adjust to your path
#OUTPUT_DIR = "cafe_filtered_orthogroups"
CAFE_COUNT_FILE = "cafe5_gene_counts.tsv"

#os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load orthogroups
df = pd.read_csv(ORTHOGROUPS_TSV, sep="\t", index_col=0)
print(f"Total orthogroups: {len(df)}")

def count_genes(cell):
    if pd.isna(cell) or str(cell).strip() == '':
        return 0
    return len(str(cell).split(', '))

count_df = df.map(count_genes)
#print(count_df)
# Standard CAFE5 filters
filtered = count_df[
    # Present in at least half your species (10 of 21)
    (count_df.gt(0).sum(axis=1) >= 10) &
    # No single species has an absurdly large count (likely annotation artifact)
    (count_df.max(axis=1) <= 100) &
    # At least some variation exists (not identical counts everywhere)
    (count_df.std(axis=1) > 0)
]
print(f"After standard filters: {len(filtered)}")

# Remove orthogroups where every present species has exactly 1 copy
variable = filtered[filtered.apply(
    lambda row: not all(v == 1 for v in row if v > 0), axis=1
)]
print(f"After removing invariant single-copy: {len(variable)}")

# Copy filtered orthogroup FASTA files to output directory
#print(f"\nCopying {len(variable)} orthogroup files to {OUTPUT_DIR}...")
#copied = 0
#missing = 0

#for og_id in variable.index:
    # OrthoFinder names files as OG0000001.fa or OG0000001.faa
#    for ext in ['.fa', '.faa', '.fasta']:
#        src = os.path.join(ORTHOGROUP_SEQS_DIR, f"{og_id}{ext}")
#        if os.path.exists(src):
#            dst = os.path.join(OUTPUT_DIR, f"{og_id}{ext}")
#            shutil.copy2(src, dst)
#            copied += 1
#            break
#    else:
#        missing += 1

#print(f"  Copied: {copied}")
#print(f"  Missing: {missing}")
#if missing > 0:
#    print(f"  Note: missing files may mean ORTHOGROUP_SEQS_DIR path needs adjusting")

# Generate CAFE5 count matrix
# CAFE5 format: Desc <tab> ID <tab> Species1 <tab> Species2 ...
print(f"\nWriting CAFE5 count matrix to {CAFE_COUNT_FILE}...")
cafe_df = variable.copy()
cafe_df.insert(0, 'Desc', '(null)')
cafe_df.index.name = 'ID'
cafe_df.reset_index(inplace=True)
cafe_df.to_csv(CAFE_COUNT_FILE, sep="\t", index=False)
print(f"  Written: {len(cafe_df)} orthogroups x {len(variable.columns)} species")

print(f"\nDone!")
#print(f"  Filtered FASTA files: {OUTPUT_DIR}/")
print(f"  CAFE5 input matrix:   {CAFE_COUNT_FILE}")
print(f"\nNext: run CAFE5 with an ultrametric tree:")
print(f"  cafe5 -i {CAFE_COUNT_FILE} -t ultrametric_tree.nwk -o cafe5_output/")