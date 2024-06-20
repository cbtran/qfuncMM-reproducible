#!/bin/bash

specs=("std" "fgn" "ar2" "anisotropic", "diag_time")

mkdir -p ./full-run/data

# If no command line arguments were provided, default to "all"
if [ $# -eq 0 ]; then
    set -- "all"
fi

# Check if the input is one of the possible specs or "all"
if [[ " ${specs[*]} " == *" $1 "* ]] || [[ "$1" == "all" ]]; then
    # If "all", loop through all possible specs
    if [[ "$1" == "all" ]]; then
        for setting in "${specs[@]}"; do
            for delta in low mid high; do
                for psi in low mid high; do
                    Rscript R_files/generate.R $delta $psi $setting 100 ${2:-1234}
                done
            done
        done
    else
        for delta in low mid high; do
            for psi in low mid high; do
                Rscript R_files/generate.R $delta $psi $1 100 ${2:-1234}
            done
        done
    fi
else
    echo "Invalid specification. Please provide one of the following: std, fgn, ar2, anisotropic. Run with no arguments to generate all specifications."
fi