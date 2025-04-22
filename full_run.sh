#!/bin/bash

# Example usage:
# ./full_run.sh std high high 1 standard
# ./full_run.sh std all 1e-2 noiseless

specs=("std" "fgn" "ar2" "anisotropic" "diag_time")
levels=("low" "mid" "high")

mkdir -p ./full-run/out

# Check if $1 is one of the possible specs
if [[ " ${specs[*]} " != *" $1 "* ]]; then
    echo "Invalid specification. Please provide one of the following: std, fgn, ar2, anisotropic."
    exit 1
fi

if [[ " ${levels[*]} " != *" $2 "* ]] && [[ "$2" != "all" ]]; then
    echo "Invalid delta. Please provide one of the following: low, mid, high, all."
    exit 1
fi

if [[ " ${levels[*]} " != *" $3 "* ]] && [[ "$2" != "all" ]]; then
    echo "Invalid psi. Please provide one of the following: low, mid, high."
    exit 1
fi

# If $2 is "all", loop through all possible levels for both delta and psi
if [[ "$2" == "all" ]]; then
    for delta in "${levels[@]}"; do
        for psi in "${levels[@]}"; do
            Rscript full-run/full-run.R $1 $delta $psi $2 $3 1 100
        done
    done
else
    Rscript full-run/full-run.R $1 $2 $3 $4 $5 1 100
fi
