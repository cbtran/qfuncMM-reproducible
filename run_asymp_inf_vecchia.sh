delta=$1
psi=$2

if [ -z "$delta" ] || [ -z "$psi" ]; then
    echo "Usage: $0 <delta> <psi>"
    exit 1
fi

nohup Rscript R_files/analysis/compute_asymp_ci.R out std stage2_vecchia vecchia $delta $psi > vec-$delta-$psi-asymp.log 2>&1 &