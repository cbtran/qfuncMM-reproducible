#!/bin/bash
# Run with ./generate_data.sh nsim seed
# nsim is optional and defaults to 100
# seed is optional and defaults to 1234

export BLAS_NUM_THREADS=1
export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

specs=("std" "fgn" "ar2" "anisotropic")

outdir=$1
mkdir -p $outdir

nsim=${2:-100}
seed=${3:-1234}

noise_level=1
for setting in "${specs[@]}"; do
    echo "Generating $setting data with noise level $noise_level"
    for delta in low mid high; do
        for psi in low mid high; do
            suffix="$setting/"
            output_dir="$outdir/$suffix"
            mkdir -p "$output_dir"
            Rscript R_files/generate.R $delta $psi $setting $noise_level $nsim $seed "$output_dir"
        done
    done
done

noise_levels=("1e-2" "1e-3" "1e-7")
setting="std"
for noise_level in "${noise_levels[@]}"; do
    echo "Generating $setting data with noise level $noise_level"
    for delta in low mid high; do
        for psi in low mid high; do
            suffix="$setting-noise-$noise_level/"
            output_dir="$outdir/$suffix"
            mkdir -p "$output_dir"
            Rscript R_files/generate.R $delta $psi $setting $noise_level $nsim $seed "$output_dir"
        done
    done
done
