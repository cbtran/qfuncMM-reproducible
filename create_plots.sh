#!/bin/bash
Rscript R_files/visualization/create_correctly_specified_plot.R
Rscript R_files/visualization/create_hcp_sim_plot.R
Rscript R_files/visualization/create_compare_vecchia_plot.R
Rscript R_files/visualization/create_misspecified_plot.R out ar2 plots
Rscript R_files/visualization/create_misspecified_plot.R out anisotropic plots
Rscript R_files/visualization/create_misspecified_plot.R out fgn plots
Rscript R_files/visualization/create_misspecified_plot.R out nonsep_matern plots