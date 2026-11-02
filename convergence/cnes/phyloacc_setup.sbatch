#!/bin/bash
#SBATCH -J phyloacc_setup
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64        # match num_procs in the config
#SBATCH --mem=128G                # whole concatenated alignment is read into memory
#SBATCH -t 48:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH --mail-user=zoe.muller@uta.edu
#SBATCH --mail-type=END,FAIL
#
# Step 1: PhyloAcc's interface reads the alignment, filters loci, builds gene
# trees (theta), splits loci into batches and writes the Snakemake workflow.
# It ends by printing a snakemake command with --dryrun.
module load biocontainers/default
module load anaconda
ml openjdk
conda activate phyloacc            # env with phyloacc, iqtree2, java (+ ASTRAL jar)
cd /anvil/scratch/x-zlm9173

phyloacc.py --config phyloacc_config.yaml
