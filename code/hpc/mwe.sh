#!/bin/bash
#PBS -j oe
#PBS -N mwe
#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=1:mem=1gb

module load tools/prod
module load R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
Rscript sim/mwe.r
