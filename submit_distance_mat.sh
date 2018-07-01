#!/bin/csh

#$ -N dist_
#$ -t 1

module load OpenMPI/.1.10.2/R/3.3.3

Rscript distance_mat.R $SGE_TASK_ID