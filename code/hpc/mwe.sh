#!/bin/bash
#PBS -j oe
#PBS -N mwe
#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=1:mem=4gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/mwe.r
