#!/bin/bash
#PBS -j oe
#PBS -N cseb.hom
#PBS -l walltime=0:45:00
#PBS -l select=1:ncpus=128:mem=512gb
#PBS -J 1-10

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/cseb.r .k=hom .cores=128 .nb=10 .b=$PBS_ARRAY_INDEX
