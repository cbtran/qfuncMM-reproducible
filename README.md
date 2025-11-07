# Reproduce qfuncMM simulation results

This repository contains R code to reproduce the simulation results in [A Mixed Model Approach for Estimating Regional Functional Connectivity from Voxel-level BOLD Signals](https://arxiv.org/abs/2211.02192).

## Generation

The script to generate simulations is `R_files/generate.R` and should be called
from the terminal. For example, the following will generate 100
simulations using a high "delta" and mid "phi" setting with correct specification.

```
Rscript R_files/generate.R high mid std 100
```

Available specifictions are `std` and misspecifications `ar2`, `fgn`, `anisotropic`.

## Install package

Our `qfuncMM` package is available on [github](https://github.com/cbtran/qfuncMM), which can be install using `devtools`:

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("cbtran/qfuncMM")
```

## Run

To fit the mixed model, run the script `full-run/full-run.R` after generating the data
and installing the `qfuncMM` package.
For example run the following after generating with the above settings:

```
Rscript full-run/full-run.R high mid std 1 100 > high-mid-std.log
```

The above pipes output from the optimization process to a log file.
Note that running 100 simuations will take quite some time.
Please see `full-run/full-run.html` for a detailed description of our simulation settings and full simulation results.
