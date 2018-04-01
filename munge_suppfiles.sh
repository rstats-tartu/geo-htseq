#!/bin/bash

#SBATCH --job-name=munge_suppfiles
#SBATCH --workdir=/gpfs/hpchome/taavi74/Projects/geo-htseq-article
#SBATCH --mail-user=tapa741@gmail.com
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --mem=36G
#SBATCH --time=08:00:00
#SBATCH --array=1-2

module purge
module load R-3.4.1

ARGS[1]=TRUE
ARGS[2]=FALSE

srun Rscript R/06_munge_suppfiles.R ${ARGS[$SLURM_ARRAY_TASK_ID]}
