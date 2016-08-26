#!/bin/bash

#PBS -l vmem=5G,walltime=23:00:00

home=/hpf/largeprojects/agoldenb/dustin/Data/TCGA/2016_01_28

module load R/3.2.3

cancer=$CANCER


cd $home/$cancer/genome-wide

Rscript --vanilla clustering-runs.R --cancer $cancer