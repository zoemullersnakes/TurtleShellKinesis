import os
import re
from Bio import SeqIO


def _base_accession(record_id, species):
    """Strip a leading 'Species_' prefix (added for Tiberius-derived CDS
    files to keep gene IDs unique across species) and any trailing numeric
    suffixes, so 'Clemmys_guttata_g2379.t1' and 'g2379' both reduce to the
    same base gene ID 'g2379'. Repeats the strip because gffread appends an
    extra '.N' to disambiguate duplicate transcript IDs (e.g. Tiberius
    restarting gene numbering per scaffold produces two different 'g1.t1'
    loci, which gffread writes out as 'g1.t1' and 'g1.t1.1')."""
    base = record_id
    prefix = f"{species}_"
    if base.startswith(prefix):
        base = base[len(prefix):]
    while True:
        new_base = re.sub(r'(\.t\d+|_t\d+|\.\d+)$', '', base)
        if new_base == base:
            break
        base = new_base
    return base


def cds_orthologs(input_cds_dir, input_orthologs_dir, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    # Step 1: load all CDS sequences into memory, keyed by species then accession.
    # We keep two indexes per species:
    #   exact_db[species][record_id]      -> exact header match (NCBI-style XP_ etc.)
    #   base_db[species][base_gene_id]    -> fallback match on bare gene ID
    #                                        (handles species-prefixed / transcript-suffixed
    #                                        headers from Tiberius-derived CDS files)
    exact_db = {}
    base_db = {}
    for f in os.listdir(input_cds_dir):
        if not f.endswith("_CDS_byprotein.fasta") or f.startswith("."):
            continue
        species = f.replace("_CDS_byprotein.fasta", "")
        exact_db[species] = {}
        base_db[species] = {}
        path = os.path.join(input_cds_dir, f)
        for rec in SeqIO.parse(path, "fasta"):
            exact_db[species][rec.id] = rec
            base = _base_accession(rec.id, species)
            if base in base_db[species]:
                print(f"  WARNING: {species} has multiple CDS records reducing to base ID "
                      f"'{base}' ({base_db[species][base].id} vs {rec.id}) - "
                      f"only the first-loaded one will be used for fallback matches")
            else:
                base_db[species][base] = rec
        print(f"Loaded {len(exact_db[species])} CDS seqs for {species}")

    # Case-insensitive species lookup: maps lowercased species name -> the
    # canonical (as-cased-in-filename) species key used in exact_db/base_db,
    # so "chrysemys_picta_bellii" in an ortholog header still finds the
    # "Chrysemys_picta_bellii" CDS file.
    species_lookup = {sp.lower(): sp for sp in exact_db}

    # Step 2: process each ortholog alignment file
    ortholog_files = [f for f in os.listdir(input_orthologs_dir)
                      if f.endswith(".fa") and not f.startswith(".")]

    for f in ortholog_files:
        og_name = f.replace(".fa", "")
        input_path = os.path.join(input_orthologs_dir, f)
        output_path = os.path.join(output_folder, f)

        recs_to_write = []
        missing = []

        for record in SeqIO.parse(input_path, "fasta"):
            parts = record.id.split("-", 1)  # split only on the FIRST dash
            if len(parts) != 2:
                print(f"  WARNING: unexpected header format: {record.id}")
                continue
            species_raw, rest = parts
            accession = rest.split("|")[0]  # e.g. "g2379|g2379_t1|contig_..." -> "g2379"

            species = species_lookup.get(species_raw.lower())
            if species is None:
                missing.append(f"no CDS file for {species_raw}")
                continue

            cds_rec = exact_db[species].get(accession) or base_db[species].get(accession)
            if cds_rec is None:
                missing.append(f"{accession} not found in {species}")
                continue

            # Get the CDS sequence and rename header to Species-Accession
            # (using the CDS file's canonical casing, not the ortholog header's)
            cds_rec = cds_rec[:]  # shallow copy so we don't mutate the shared cached record
            cds_rec.id = f"{species}-{accession}"
            cds_rec.description = ""
            recs_to_write.append(cds_rec)

        # NOTE: this block must stay INSIDE the "for f in ortholog_files" loop,
        # so each orthogroup gets its own output file.
        if recs_to_write:
            SeqIO.write(recs_to_write, output_path, "fasta")
            print(f"  {og_name}: wrote {len(recs_to_write)} sequences")
        else:
            print(f"  {og_name}: WARNING - no sequences written")

        if missing:
            for m in missing[:3]:  # show first 3 warnings only
                print(f"    MISSING: {m}")


# Run it
cds_orthologs(
    input_cds_dir="/media/user/Draco/turtles/protein_files",
    input_orthologs_dir="/media/user/Draco/turtles/protein_files/proteins/primary_transcripts/OrthoFinder/Results_Sep08/Single_Copy_Orthologue_Sequences",
    output_folder="/media/user/Draco/turtles/protein_files/singlecopy_nuc"
)
