#!/bin/bash
#PBS -j oe
#PBS -N depr
#PBS -l walltime=0:30:00
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -J 1-10000

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/depr.r \
  .cores=1 .debug=0 .nb=10000 .b=$PBS_ARRAY_INDEX

# total runs = 216090  <- 21 * 5 * 6 * 7 * 7
# batch runs = 22      <- 216090 / 10000
# batch time = 0:30:00 <- 22 * 0:30 / 1 cpu = 0:11:00 (+)
# batch mem  = 6GB     <- 1.5 GB (+)
