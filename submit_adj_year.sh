#!/bin/csh

#$ -N adj_
#$ -t 1-6

module load bio/R/3.3.1-gcc

Rscript adj_year.R $SGE_TASK_ID