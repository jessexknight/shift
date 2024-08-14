#!/bin/bash
#PBS -j oe
#PBS -N val
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=64:mem=128gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR

for v in {1..23}; do
  Rscript sim/val.r v=$v n=1000 n.seed=128 .cores=64 .mem=128
done
