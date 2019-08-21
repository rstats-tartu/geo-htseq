#!/bin/sh
#SBATCH --job-name=scopus
#SBATCH --partition=main
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000
#SBATCH --output=slurm_%j.out

module load singularity \
	&& cd $HOME/Projects/geo-rnaseq \
	&& singularity exec $HOME/simg/singularity-tidyverse.simg Rscript R/download_scopus_citations.R

