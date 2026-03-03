#!/bin/bash
#PBS -j oe
#PBS -N lhs.RRo
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-10000

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/mass.emu.r \
  .debug=0 .k=RRo .cores=1 .nb=10000 .b=$PBS_ARRAY_INDEX
