#!/bin/bash
#PBS -j oe
#PBS -N cseb
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=128:mem=256gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/cseb.r .cores=128
