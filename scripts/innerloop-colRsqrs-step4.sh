#!/bin/bash

#PBS -l nodes=1:ppn=30,vmem=5G,walltime=99:00:00

home=/hpf/largeprojects/agoldenb/dustin/Data/TCGA/2016_01_28

module load R/3.2.3

cancer=$CANCER
k=$K


cd $home/$cancer/genome-wide

Rscript --vanilla calling-rsqr-generation-functions.r --cancer $cancer -k $k -p T -t two -l T --cores 30 -s four

cd $home
qsub innerloop-colRsqrs-step4.sh -v CANCER=$cancer,K=$k
