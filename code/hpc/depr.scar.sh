#!/bin/bash
#PBS -j oe
#PBS -N depr.scar
#PBS -l walltime=1:30:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-1000

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/app/depr.r \
  .k=scar .cores=1 .debug=0 .nb=1000 .b=$PBS_ARRAY_INDEX

# total runs = 61740   <- 21 * 5 * 6 * 7 * 14
# batch runs = 62      <- 61740 / 1000
# batch time = 1:30:00 <- 62 * 0:30 / 1 cpu = 0:31:00 (+)
# batch mem  = 4GB     <- 1.5 GB (+)
