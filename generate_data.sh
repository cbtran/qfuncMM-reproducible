#!/bin/bash
# Run with ./generate_data.sh data_dir nsim seed
# nsim is optional and defaults to 100
# seed is optional and defaults to 1234
# This script generates only the correctly specified setting (std with noise_level=1)

export BLAS_NUM_THREADS=1
export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

outdir=$1
mkdir -p $outdir

nsim=${2:-100}
seed=${3:-1234}

# Generate only the correctly specified setting
setting="std"
noise_level=1
echo "Generating $setting data with noise level $noise_level"
for delta in low mid high; do
    for psi in low mid high; do
        suffix="$setting/"
        output_dir="$outdir/$suffix"
        mkdir -p "$output_dir"
        Rscript R_files/simulation/generate.R $delta $psi $setting $noise_level $nsim $seed "$output_dir"
    done
done
