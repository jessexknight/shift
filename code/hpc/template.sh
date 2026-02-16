#!/bin/bash
#PBS -j oe
#PBS -N <<NAME>>
#PBS -l walltime=<<HMS>>
#PBS -l select=1:ncpus=<<CPU>>:mem=<<MEM>>gb
#PBS -J 1-<<NB>>

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript <<SCRIPT>> \
  <<FLAGS>> .cores=<<CPU>> .nb=<<NB>> .b=$PBS_ARRAY_INDEX
