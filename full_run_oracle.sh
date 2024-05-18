#!/bin/bash

levels=("low" "mid" "high")

mkdir -p ./full-run/out

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
            Rscript R_files/full-run/full-run-oracle.R $delta $psi 1 100
        done
    done
else
    Rscript R_files/full-run/full-run-oracle.R $2 $3 1 100
fi