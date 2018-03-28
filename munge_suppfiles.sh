#!/bin/bash

#The name of the job is test_job
#SBATCH -J munge_suppfiles

#The job requires N compute nodes
#SBATCH --ntasks=1

# Allocate some good amount of memory
#SBATCH --mem=45G

# Send email when finished
#SBATCH --mail-type=END
#SBATCH --mail-user=tapa741@gmail.com

#SBATCH --array=1-2

#These commands are run on one of the nodes allocated to the job (batch node)
module purge
module load R-3.4.1

cd /gpfs/hpchome/taavi74/Projects/geo-htseq-article

ARGS=(TRUE FALSE)

srun Rscript R/06_munge_suppfiles.R ${ARGS[$SLURM_ARRAY_TASK_ID]}

