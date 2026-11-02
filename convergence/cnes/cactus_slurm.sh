#!/bin/bash

#SBATCH -J partial_turtle_cactus
#SBATCH -N 1
#SBATCH -p shared
#SBATCH -n 64
#SBATCH -t 72:00:00
#SBATCH --mail-user=zoe.muller@uta.edu
#SBATCH --mail-type=all

module --force purge
ml biocontainers cactus

cactus /anvil/scratch/x-zlm9173/jobStore partial.txt partial.hal --batchSystem slurm --consCores 64 --doubleMem true --batchLogsDir /anvil/scratch/x-zlm9173/cactus_jobs --slurmTime 72:00:00