#!/bin/bash
# Run with ./generate_data.sh std noise_level seed
# seed is optional and defaults to 1234

specs=("std" "fgn" "ar2" "anisotropic")

mkdir -p ./full-run/data

# If no command line arguments were provided, default to "all"
if [ $# -eq 0 ]; then
    set -- "all"
fi

# Set variables from arguments with defaults
spec=$1
noise_level=${2:-1}
seed=${3:-1234}

# Check if the input is one of the possible specs or "all"
if [[ " ${specs[*]} " == *" $spec "* ]] || [[ "$spec" == "all" ]]; then
    # If "all", loop through all possible specs
    if [[ "$spec" == "all" ]]; then
        for setting in "${specs[@]}"; do
            for delta in low mid high; do
                for psi in low mid high; do
                    Rscript R_files/generate.R $delta $psi $setting $noise_level 100 $seed
                done
            done
        done
    else
        for delta in low mid high; do
            for psi in low mid high; do
                Rscript R_files/generate.R $delta $psi $spec $noise_level 100 $seed
            done
        done
    fi
else
    echo "Invalid specification. Please provide one of the following: std, fgn, ar2, anisotropic. Run with no arguments to generate all specifications."
fi