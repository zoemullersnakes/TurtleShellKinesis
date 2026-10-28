# Need to annotate Graptemys geographica, Deirochelys reticularia, Glyptemys, Clemmys, Emydoidea, Actinemys, Pangshura

# Run: conda activate tiberius 

import subprocess
import os

""" Genomes that need annotations:
- Graptemys geographica
- Deirochelys reticularia
- Glyptemys insculpta
- Clemmys guttata
- Emydoidea blandingii
- Actinemys pallida
- Pangshura tecta
"""
genomes_to_annotate = ['/media/user/Rhinoderma/turtes/new_genomes/graptemys_geographica.fna',
                       '/media/user/Rhinoderma/turtes/new_genomes/deirochelys_reticularia.fna',
                       '/media/user/Rhinoderma/turtes/new_genomes/glyptemys_insculpta.fna',
                       '/media/user/Rhinoderma/turtes/new_genomes/clemmys_guttata.fna',
                       '/media/user/Rhinoderma/turtes/new_genomes/emydoidea_blandingii.fna',
                       '/media/user/Rhinoderma/turtes/new_genomes/actinemys_pallida.fna',
                       '/media/user/Rhinoderma/turtes/new_genomes/pangshura_tecta.fna']
output_path = "/media/user/Draco/turtles/protein_files/"

# Ensure output directory exists
os.makedirs(output_path, exist_ok=True)

for genome in genomes_to_annotate:
    print(f"Starting tiberius on {genome}")

    # Extract base name using os.path.basename
    base = os.path.splitext(os.path.basename(genome))[0]

    # Construct the command as a list of strings for subprocess.run
    command = [
        "python",
        "tiberius.py",
        "--genome", genome,
        "--model_cfg", "vertebrates",
        "--out", os.path.join(output_path, f"{base}.gff3"),
        "--protseq", os.path.join(output_path, f"{base}.pep"),
        "--codingseq", os.path.join(output_path, f"{base}.cds")
    ]

    try:
        # Execute the command using subprocess.run
        # You might need to adjust the path to 'tiberius.py' if it's not in the current directory
        # If 'tiberius.py' needs to be run from a specific conda environment, ensure it's activated before running this script
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running tiberius on {genome}:\n{e.stdout}\n{e.stderr}")
    except FileNotFoundError:
        print("Error: 'python' or 'tiberius.py' command not found. Ensure Python and the script are in your PATH.")

    print(f"Finished tiberius on {genome}")
