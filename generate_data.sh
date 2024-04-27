#!/bin/bash

settings=("std" "fgn" "ar2" "anisotropic")

mkdir -p ./full-run/data

# Check if the input is one of the possible settings or "all"
if [[ " ${settings[*]} " == *" $1 "* ]] || [[ "$1" == "all" ]]; then
    # If "all", loop through all possible settings
    if [[ "$1" == "all" ]]; then
        for setting in "${settings[@]}"; do
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
    echo "Invalid specification. Please provide one of the following: std, fgn, ar2, anisotropic, all."
fi