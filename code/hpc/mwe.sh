#!/bin/bash
#PBS -lwalltime=00:01:00
#PBS -lselect=1:ncpus=1:mem=1gb

module load tools/prod
module load R/4.1.2-foss-2021b

cd $PBS_O_WORKDIR
Rscript sim/mwe.r
