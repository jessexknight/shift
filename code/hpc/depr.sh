#!/bin/bash
#PBS -j oe
#PBS -N depr
#PBS -l walltime=0:30:00
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -J 1-100

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/depr.r \
  .cores=16 .debug=0 .k=hetc .nb=100 .b=$PBS_ARRAY_INDEX
