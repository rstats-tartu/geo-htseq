#!/bin/bash

#SBATCH --job-name=munge_suppfiles
#SBATCH --workdir=/gpfs/hpchome/taavi74/Projects/geo-htseq-article
#SBATCH -p main
#SBATCH --mail-user=tapa741@gmail.com
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --array=1-100
#SBATCH -o output/slurm_log/slurm-%A_%a.out

module purge
module load R-3.4.1

srun Rscript R/06_munge_suppfiles_slurm.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX
