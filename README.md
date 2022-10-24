# qfuncMM-reproduce

This repository contains R codes to reproduce the simulation results in our paper [Quantifying brain functional connectivity with noisy voxel-level signals](https://arxiv.org/) paper.

## Overview

All R scripts are included in folder `R_files`. The bash file to reproduce all simulation is included in folder `bash`. Simulation data will be saved in folder `simulation-data`, and simulation outputs will be saved in folder `simulation-output`. You can modify `R_files/simulation.R` to change the location of data and output.

## Install package

Our `qfuncMM` package is available on [github](https://github.com/cbtran/qfuncMM), which can be install using `devtools`:

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("cbtran/qfuncMM")
```

To install dependencies for simulation, type in terminal console while inside `qfuncMM-reproducible` folder

```
Rscript ./R_files/install_dependencies.R
```


# Replicate 1 simulation setting

Each script in `simulation_setting` contains information for 1 simulation setting. You can specify signal strength, intra-correlation, number of voxels in each region, side length of 3D lattice, number of timepoints, number of simulation, and number of B-spline basis as mentioned in the manuscript. 
The order of command line input is:
1. k_eta
2. phi_gamma
3. number of voxels
4. side length of 3D lattice. For example, a 3D lattice with side length 5 has total 5*5*5 voxels.
5. number of timepoints
6. number of replications
7. number of bspline basis

For example, to run 10 simulations with weak signal `k_eta = 0.5`, weak intra-correlation `phi_gamma = 1`, 50 voxels `L=50`, side length 7, 60 timepoints `M=60`, and 45 B-splines basis, run the following in terminal console while inside `qfuncMM-reproducible` folder
 
```
Rscript ./R_files/simulation.R 0.5 1 50 7 60 10 45
```

# Replicate all simulation setting

To replicate all simulation, you can use bash file `bash/reproduce_all.sh` as follow

```
bash ./bash/reproduce_all.sh
```

Note that this may take a very long time since `L=50, M=60` and number of simulation is 100.
This will produce a `Rplots.pdf` that contains all plots from the manuscript and a `results_all.csv` file that contains the table of numerical results.

