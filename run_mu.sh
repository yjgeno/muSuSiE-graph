#!/bin/bash 

trap break INT
for L_max in `seq 10 $max`
do
  echo "L_max: $L_max"
  Rscript GRNs.R $L_max
  python ./experiment/results/mu_analysis.py -L $L_max
  echo -e "run $L_max completed\n"
done
trap - INT
