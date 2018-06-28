#!/bin/csh

#$ -N adj_
#$ -t 1-6

module load bio/R/3.3.1-gcc
. /usr/share/modules/init/bash

Rscript adj_year.R $SGE_TASK_ID