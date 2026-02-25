#!/bin/bash
#PBS -j oe
#PBS -N RRo.rev.base
#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-100

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/mass.r \
  .debug=0 .k=RRo.rev.base .cores=1 .nb=100 .b=$PBS_ARRAY_INDEX
