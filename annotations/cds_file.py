import re
import os
import subprocess
from Bio import SeqIO

# Configuration - adjust these paths
GFF_DIR = "/media/user/Draco/turtles/protein_files"        # folder with Species.gff/.gff3 files
GENOME_DIR = "/media/user/Rhinoderma/turtes/new_genomes"   # folder with Species_genome.fna files
OUTPUT_DIR = "/media/user/Draco/turtles/protein_files"     # where to write Species_CDS_byprotein.fasta

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Map species name to genome filename - adjust to match your actual filenames
# Key = species name (matching GFF filename without extension)
# Value = genome fasta filename
genome_files = {
	"Mauremys_reevesii":  "Mauremys_reevesii.fa",
	"Emys_orbicularis": "Emys_orbicularis.fa",
	"Chrysemys_picta_bellii": "Chrysemys_picta_bellii.fa",
	"glyptemys_insculpta":  "glyptemys_insculpta.fna",
	"actinemys_pallida":  "Actinemys_pallida.fna",
	"deirochelys_reticularia":  "deirochelys_reticularia.fna",
	"clemmys_guttata":  "Clemmys_guttata.fna",
  "Mauremys_mutica":  "Mauremys_mutica.fa",
  "Malaclemys_terrapin_pileata":  "Malaclemys_terrapin_pileata.fa",
  "Trachemys_scripta_elegans":  "Trachemys_scripta_elegans.fa",
  "Terrapene_triunguis":  "GCF_002925995.2_T_m_triunguis-2.0_genomic.fna",
  "graptemys_geographica":  "graptemys_geographica.fna",
  "pangshura_tecta":  "pangshura_tecta.fna",
  "Cuotri": "Cuora_V1.full.softmasked.fa",
  "emydoidea_blandingii": "emydoidea_blandingii.fna"
}

def find_gff(species):
    """Look for Species.gff3 first, then Species.gff (Tiberius default is .gff3,
    but can also be told to write .gff)."""
    for ext in (".gff3", ".gff"):
        candidate = os.path.join(GFF_DIR, f"{species}{ext}")
        if os.path.exists(candidate):
            return candidate
    return None


def build_gene_symbol_map(gff_path):
    """NCBI/RefSeq-style annotations: CDS lines carry a gene=<symbol> attribute
    directly (e.g. '...;gene=DOCK9;protein_id=XP_...'), in addition to
    Parent=rna-XM_.... Map transcript ID (the Parent, matching gffread's
    header) -> gene symbol. Orthogroup headers are keyed by bare gene symbol
    (e.g. 'Chrysemys_picta_bellii-DOCK9'), so this is what CDS output needs
    to match on."""
    gene_map = {}
    with open(gff_path) as f:
        for line in f:
            if "\tCDS\t" not in line:
                continue
            attrs = line.strip().split("\t")[8]
            parent_match = re.search(r'Parent=([^;]+)', attrs)
            gene_match   = re.search(r'(?:^|;)gene=([^;]+)', attrs)
            if parent_match and gene_match:
                gene_map[parent_match.group(1)] = gene_match.group(1)
    return gene_map


# Process each species
for species, genome_fname in genome_files.items():
    gff_path    = find_gff(species)
    genome_path = os.path.join(GENOME_DIR, genome_fname)
    raw_cds     = os.path.join(OUTPUT_DIR, f"{species}_CDS_raw.fasta")
    out_cds     = os.path.join(OUTPUT_DIR, f"{species}_CDS_byprotein.fasta")

    print(f"\n{'='*50}")
    print(f"Processing: {species}")

    # Check input files exist
    if gff_path is None:
        print(f"  SKIP: GFF not found: {os.path.join(GFF_DIR, species)}.gff3/.gff")
        continue
    if not os.path.exists(genome_path):
        print(f"  SKIP: genome not found: {genome_path}")
        continue

    print(f"  Using GFF: {gff_path}")

    # Step 1: run gffread to extract raw CDS
    print(f"  Running gffread...")
    result = subprocess.run(
        ["gffread", gff_path, "-g", genome_path, "-x", raw_cds],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  ERROR: gffread failed: {result.stderr}")
        continue
    if not os.path.exists(raw_cds) or os.path.getsize(raw_cds) == 0:
        print(f"  ERROR: gffread produced empty output")
        continue

    # Step 2: try to build a transcript -> gene symbol map (NCBI/RefSeq-style
    # annotations). Tiberius (and other ab initio predictors like AUGUSTUS)
    # don't carry a gene= attribute on CDS lines, so this comes back empty
    # and we fall back to gffread's own transcript IDs.
    print(f"  Checking GFF ID style...")
    gene_map = build_gene_symbol_map(gff_path)

    if gene_map:
        # ---- NCBI-style path: one sequence per gene, keyed by gene symbol ----
        print(f"  Detected NCBI-style GFF: mapped {len(gene_map)} transcripts -> gene symbols")

        # A gene can have multiple transcripts/CDS but the orthogroup headers
        # carry no isoform info, so keep only the longest CDS per gene symbol
        # (the usual "primary transcript" convention).
        best_by_gene = {}  # gene_symbol -> SeqRecord
        unmapped = 0
        for rec in SeqIO.parse(raw_cds, "fasta"):
            gene = gene_map.get(rec.id)
            if gene is None:
                unmapped += 1
                continue
            current_best = best_by_gene.get(gene)
            if current_best is None or len(rec.seq) > len(current_best.seq):
                best_by_gene[gene] = rec

        recs = []
        for gene, rec in best_by_gene.items():
            rec.id = gene
            rec.description = ""
            recs.append(rec)

        SeqIO.write(recs, out_cds, "fasta")
        print(f"  Wrote {len(recs)} sequences (longest isoform per gene), {unmapped} transcripts had no gene= match")

    else:
        # ---- Tiberius / ab-initio path ----
        # No gene symbol exists; just take gffread's transcript IDs as-is
        # (e.g. g1.t1, g2.t1, ...). These IDs are NOT unique across species
        # (every Tiberius run starts numbering at g1), so prefix with species
        # name to keep downstream cross-species files sane.
        print(f"  No gene= attribute found on CDS lines - treating as ab initio (Tiberius-style) GFF")
        recs = []
        for rec in SeqIO.parse(raw_cds, "fasta"):
            rec.id = f"{species}_{rec.id}"
            rec.description = ""
            recs.append(rec)

        SeqIO.write(recs, out_cds, "fasta")
        print(f"  Wrote {len(recs)} sequences (prefixed with species name)")

    # Clean up raw CDS file
    os.remove(raw_cds)

print(f"\nDone! CDS files written to {OUTPUT_DIR}")