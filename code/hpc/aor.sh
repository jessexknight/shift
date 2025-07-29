#!/bin/bash
#PBS -j oe
#PBS -N aor.o6
#PBS -l walltime=0:45:00
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -J 1-100

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/aor.r .k=o6 .cores=16 .nb=100 .b=$PBS_ARRAY_INDEX
