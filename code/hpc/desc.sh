#!/bin/bash
#PBS -j oe
#PBS -N desc
#PBS -l walltime=0:30:00
#PBS -l select=1:ncpus=64:mem=512gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/desc/dep.r .cores=64
