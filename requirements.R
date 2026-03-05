installed <- rownames(installed.packages())

if (!"fields" %in% installed) {
  install.packages("fields", repos = "https://cloud.r-project.org")
}

if (!"multiwave" %in% installed) {
  install.packages("multiwave", repos = "https://cloud.r-project.org")
}

if (!"waveslim" %in% installed) {
  install.packages("waveslim", repos = "https://cloud.r-project.org")
}

if (!"jsonlite" %in% installed) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}

if (!"nleqslv" %in% installed) {
  install.packages("nleqslv", repos = "https://cloud.r-project.org")
}

if (!"remotes" %in% installed) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

if (!"GpGpQFuncMM" %in% installed) {
  # Note: This package requires gfortran/Fortran compiler to be installed
  # On macOS with Homebrew: brew install gfortran
  # On Ubuntu/Debian: sudo apt-get install gfortran
  # On Windows: Install Rtools which includes gfortran
  # A custom Makevars file at ~/.R/Makevars may be needed for proper linking
  remotes::install_github("roobnloo/GpGp-qfuncMM")
}

if (!"qfuncMM" %in% installed) {
  remotes::install_github("cbtran/qfuncMM")
}
