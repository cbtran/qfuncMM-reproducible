#!/bin/bash
signal_list=("0.5" "1" "1.5")
intra_list=("1" "1/2" "1/4")

for signal in ${signal_list[@]}
do
  for intra in ${intra_list[@]}
  do
    Rscript ./R_files/simulation.R $signal $intra 50 7 60 100 45
  done
done
