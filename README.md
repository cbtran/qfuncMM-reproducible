# Reproduce qfuncMM simulation results

This repository contains R code to reproduce the simulation results in our
upcoming paper.  An older draft is available [here](https://arxiv.org/abs/2211.02192).

## Install packages

The file `requirements.R` lists the set of packages that implement our method as well as auxiliary packages used for analysis.
All dependencies may be installed by running in the terminal

```
> Rscript requirements.R
```


## Data generation

The script to generate simulation data is `generate_data.sh`. For example, the following will generate 100 simulated datasets in each setting and save them to the `data/` directory. The total size for all the data is around `400MB`.

```
./generate_data.sh data 100
```

## Run stage 1

First, stage 1 is run on all datasets to estimate the intra-regional parameters for each simulation setting. This can be run with `run_stage1.sh`.
For example, the following will run the noisy stage 1 model under the correctly specified setting, reading the simulation data from `data/` and outputting the results to `out/`.

```
./run_stage1.sh std noisy data out
```

## Run stage 2

After running stage 2, stage 2 can be run with the `run_stage2_*.sh` scripts. To run the ReML model, `run_stage2_reml.sh` is used whereas `run_stage2_vecchia.sh` is used for the Vecchia's approximation.
For instance, the following will run Vecchia's approximation of stage 2 under the previously mentioned stage 1 setting and save the results in the corresponding subfolder in `out/`.

```
./run_stage2_vecchia.sh std noisy data out
```
