#!/bin/bash
#PBS -j oe
#PBS -N depr
#PBS -l walltime=0:20:00
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -J 1-1000

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/depr.r \
  .cores=16 .debug=0 .nb=1000 .b=$PBS_ARRAY_INDEX

# total runs = 216090  <- 21 * 5 * 6 * 7 * 7
# batch runs = 216     <- 216090 / 1000
# batch time = 0:20:00 <- 216 * 1:00 (+) / 16 cpu = 0:13:50 (+)
# batch mem  = 64      <- 3 GB (+) * 16 = 48 GB (+)
