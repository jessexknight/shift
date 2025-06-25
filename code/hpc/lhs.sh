#!/bin/bash
#PBS -j oe
#PBS -N lhs
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -J 1-9

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/lhs.r .cores=128 .nb=9 .b=$PBS_ARRAY_INDEX
