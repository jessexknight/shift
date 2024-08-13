#!/bin/bash
#PBS -lwalltime=00:10:00
#PBS -lselect=1:ncpus=64:mem=8gb

module load tools/prod
module load R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
Rscript sim/val.r n=1000 n.seed=128 .cores=128
