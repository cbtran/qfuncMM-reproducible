#!/bin/bash

settings=("std" "fgn" "ar2" "anisotropic")

mkdir -p ./full-run/out

# Check if the input is one of the possible settings or "all"
if [[ " ${settings[*]} " == *" $1 "* ]] || [[ "$1" == "all" ]]; then
    # If "all", loop through all possible settings
    if [[ "$1" == "all" ]]; then
        for setting in "${settings[@]}"; do
            for delta in low mid high; do
                for psi in low mid high; do
                    Rscript R_files/full-run/full-run.R $delta $psi $setting 1 100
                done
            done
        done
    else
        for delta in low mid high; do
            for psi in low mid high; do
                Rscript R_files/full-run/full-run.R $delta $psi $setting 1 100
            done
        done
    fi
else
    echo "Invalid specification. Please provide one of the following: std, fgn, ar2, anisotropic, all."
fi