#!/bin/csh

#$ -N dist_
#$ -t 1

module load bio/R/3.3.1-gcc
. /usr/share/modules/init/bash

Rscript distance_mat.R $SGE_TASK_ID
