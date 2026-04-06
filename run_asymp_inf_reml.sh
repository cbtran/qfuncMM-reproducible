echo Creating CIs for ReML...
for delta in low mid high; do
    for psi in low mid high; do
        Rscript R_files/analysis/compute_asymp_ci.R out std stage2_reml reml $delta $psi
    done
done
echo Done