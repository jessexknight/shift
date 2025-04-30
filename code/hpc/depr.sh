#!/bin/bash
#PBS -j oe
#PBS -N depr
#PBS -l walltime=1:30:00
#PBS -l select=1:ncpus=128:mem=256gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/depr.r .cores=128
