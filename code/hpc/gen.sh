#!/bin/bash
# inputs -----------------------------------------
   HPC='hpc/mass.sh'
SCRIPT='sim/app/mass.r'
  NAME='RR2.ad.base'
 FLAGS='.debug=0 .k=RR2.ad.base'
    NR=500    # num runs (total)
   TPR=30     # time per run (sec)
   MPR=2      # mem per run (GB)
    NB=100    # num batches
   CPU=1      # num cpu (parallel)
# calcs ------------------------------------------
HMS=`date -d@$((NR*TPR/NB/CPU*2)) -u +%H:%M:%S`
MEM=$((CPU*MPR*2))
# output -----------------------------------------
cp hpc/template.sh $HPC
for i in SCRIPT NAME FLAGS NB CPU HMS MEM; do
  sed -i "s#<<$i>>#${!i}#g" $HPC
done
