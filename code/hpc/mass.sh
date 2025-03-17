#!/bin/bash
#PBS -j oe
#PBS -N mass
#PBS -l walltime=02:30:00
#PBS -l select=1:ncpus=64:mem=512gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR

/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/mass.dep.haz.r .cores=64
