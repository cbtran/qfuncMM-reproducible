#!/bin/bash
# signal_list=("0.5" "1" "1.5")
# intra_list=("1" "1/4")
# L=50 # number of voxels
# M=60 # number of timepoints
# sidelength=7 # sidelength of 3D lattice
# num_sim=100 # number of simulation
# K=45 # number of bspline basis

signal_list=("weak" "med" "strong")
intra_list=("weak" "strong")
for signal in ${signal_list[@]}
do
  for intra in ${intra_list[@]}
  do
    echo "starting signal = $signal, intra = $intra"
    Rscript ./R_files/simulation.R $signal $intra 100
    echo "finished signal = $signal, intra = $intra"
  done
done
# Rscript ./R_files/summary_all.R $M $num_sim
