export BLAS_NUM_THREADS=1
export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

# Usage: ./run_stage2_reml.sh <data_spec> <cov_setting> <data_dir> <out_dir> <delta> <psi> [oracle]
# <data_spec> is one of "std", "fgn", "ar2", "anisotropic", "nonsep", "nonsep_matern", "std-noise-1e-2", "std-noise-1e-3", "std-noise-1e-7",
# <cov_setting> is one of "noisy", "noiseless", "auto"
# <data_dir> is the directory containing the data files
# <out_dir> is the output directory
# <delta> is the delta setting ("high", "mid", "low")
# <psi> is the psi setting ("high", "mid", "low")
# [oracle] is an optional argument (TRUE or FALSE)

data_spec=$1
cov_setting=$2
data_dir=$3
out_dir=$4
delta=$5
psi=$6
oracle=${7:-FALSE} # Default to FALSE if not provided

if [ -z "$data_spec" ] || [ -z "$cov_setting" ] || [ -z "$data_dir" ] || [ -z "$out_dir" ] || [ -z "$delta" ] || [ -z "$psi" ]; then
  echo "Usage: $0 <data_spec> <cov_setting> <data_dir> <out_dir> <delta> <psi> [oracle]"
  exit 1
fi

Rscript R_files/simulation/run_stage2.R "$data_spec" "$cov_setting" "$data_dir" "$out_dir" FALSE 100 "$delta" "$psi" "$oracle"