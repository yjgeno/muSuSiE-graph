#!/bin/bash 

trap break INT
for L_max in `seq 20 $max`
do
  echo "L_max: $L_max"
  C:/Users/yjyang027/Documents/R-4.0.4/bin/Rscript GRNs.R $L_max
  python ./experiment/results/mu_analysis.py -L $L_max --file sum_unfilter
  echo -e "run $L_max completed\n"
done
trap - INT
