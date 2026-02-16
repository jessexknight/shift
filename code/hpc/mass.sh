#!/bin/bash
#PBS -j oe
#PBS -N mass
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-100

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/mass.r \
  .debug=0 .k=RRo.rev.eR2 .cores=1 .nb=100 .b=$PBS_ARRAY_INDEX
