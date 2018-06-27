#!/bin/csh

#$ -N adj_
#$ -t 1

module load bio/R/3.3.1-gcc

Rscript adj_2017.R $SGE_TASK_ID