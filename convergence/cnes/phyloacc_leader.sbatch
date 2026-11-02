#!/bin/bash
#SBATCH -J phyloacc_leader
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 96:00:00               # use shared's maximum; resubmit to resume
#SBATCH -o logs/%x_%j.out
#SBATCH --mail-user=zoe.muller@uta.edu
#SBATCH --mail-type=END,FAIL
#
# Step 2: Snakemake submits one Slurm job per batch and waits for them.
# Run only after the dry run (below) looks right. If this job hits its time
# limit, resubmit it: Snakemake skips batches that already finished.
module load biocontainers/default
module load anaconda
source ~/.bashrc
conda activate phyloacc
J=/anvil/scratch/x-zlm9173/phyloacc_out/phyloacc-job-files/snakemake

snakemake -p -s /anvil/scratch/x-zlm9173/phyloacc_out_adaptive/phyloacc-job-files/snakemake/run_phyloacc.smk \
    --configfile /anvil/scratch/x-zlm9173/phyloacc_out_adaptive/phyloacc-job-files/snakemake/phyloacc-config.yaml \
    --profile /anvil/scratch/x-zlm9173/phyloacc_out_adaptive/phyloacc-job-files/snakemake/profiles/slurm_profile \
    --rerun-incomplete --keep-going