#!/bin/bash
# Run with ./generate_additional_data.sh data_dir nsim seed
# nsim is optional and defaults to 100
# seed is optional and defaults to 1234
# This script generates additional settings beyond the correctly specified one

export BLAS_NUM_THREADS=1
export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

outdir=$1
mkdir -p $outdir

nsim=${2:-100}
seed=${3:-1234}

# Generate data for other covariance settings (fgn, ar2, anisotropic)
specs=("fgn" "ar2" "anisotropic")
noise_level=1
for setting in "${specs[@]}"; do
    echo "Generating $setting data with noise level $noise_level"
    for delta in low mid high; do
        for psi in low mid high; do
            suffix="$setting/"
            output_dir="$outdir/$suffix"
            mkdir -p "$output_dir"
            Rscript R_files/simulation/generate.R $delta $psi $setting $noise_level $nsim $seed "$output_dir"
        done
    done
done

# Generate data for different noise levels with std setting
noise_levels=("1e-2" "1e-3" "1e-7")
setting="std"
for noise_level in "${noise_levels[@]}"; do
    echo "Generating $setting data with noise level $noise_level"
    for delta in low mid high; do
        for psi in low mid high; do
            suffix="$setting-noise-$noise_level/"
            output_dir="$outdir/$suffix"
            mkdir -p "$output_dir"
            Rscript R_files/simulation/generate.R $delta $psi $setting $noise_level $nsim $seed "$output_dir"
        done
    done
done
