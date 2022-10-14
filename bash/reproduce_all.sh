#!/bin/bash
signal_list=("0.5" "1" "1.5")
intra_list=("1" "1/4")
L=50 # number of voxels
M=60 # number of timepoints
sidelength=7 # sidelength of 3D lattice
num_sim=100 # number of simulation
K=45 # number of bspline basis
for signal in ${signal_list[@]}
do
  for intra in ${intra_list[@]}
  do
    Rscript ./R_files/simulation.R $signal $intra $L $sidelength $M $num_sim $K
  done
done
Rscript ./R_files/summary_all.R $M $num_sim
