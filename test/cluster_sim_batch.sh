#!/bin/bash

#SBATCH --partition=lau
#SBATCH --array=1-3

module purge
module load R

Rscript /projects/lau_projects/phylodynamics_hannah/Within_Host_Viral_Load/batch_analysis_sim.R /projects/lau_projects/phylodynamics_hannah/Within_Host_Viral_Load/batch_parameters.csv $SLURM_ARRAY_TASK_ID