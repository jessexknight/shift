#!/bin/bash
#PBS -j oe
#PBS -N mass
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=64:mem=512gb

module add tools/prod
module add R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR

/bin/time -f "`cat hpc/time.fmt`" Rscript sim/mass/dep.haz.r .cores=64 RR.step=0.2 Ro=.03 Rx=.03
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/mass/dep.haz.r .cores=64 RR.step=0.2 Ro=.01 Rx=.03
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/mass/dep.haz.r .cores=64 RR.step=0.2 Ro=.03 Rx=.01
/bin/time -f "`cat hpc/time.fmt`" Rscript sim/mass/dep.haz.r .cores=64 RR.step=0.2 Ro=.01 Rx=.01
