suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(ggridges))

########## Read data function ########################
read_results <- function(result_folder) {
  rho_pearson <- suppressMessages(
    read_csv(file.path(result_folder, "rho_pearson.csv"),
    show_col_types = FALSE))

  rho_inter <- suppressMessages(
    read_csv(file.path(result_folder, "rho_inter.csv"),
    show_col_types = FALSE))

  rho_true <- suppressMessages(
    read_csv(file.path(result_folder, "rho_true.csv"),
    show_col_types = FALSE))

  num_sim <- nrow(rho_pearson)
  rhotype <- factor(
    rep(c("REML", "Intra", "Vanilla"), each = num_sim),
    levels = c("REML", "Intra", "Vanilla"))

  pair.12 <- data.frame(
    type = rhotype,
    rho = c(rho_inter$rho_12, rho_pearson$cor_intra_12, rho_pearson$cor_vanilla_12))

  pair.13 <- data.frame(
    type = rhotype,
    rho = c(rho_inter$rho_13, rho_pearson$cor_intra_13, rho_pearson$cor_vanilla_13))

  pair.23 <- data.frame(
    type = rhotype,
    rho = c(rho_inter$rho_23, rho_pearson$cor_intra_23, rho_pearson$cor_vanilla_23))

  pair.all <- cbind(rbind(pair.12, pair.13, pair.23),
                    inter=factor(c(rep(as.character(rho_true), each=3*num_sim)),
                                 levels=as.character(rho_true)))
  list(pair_12=pair.12, pair_13=pair.13, pair_23=pair.23, pair_all=pair.all)
}

########## Plot function ########################
box_plot_3_region <- function(result_list, true_rho) {
  
  pair_12 = result_list$pair_12
  pair_13 = result_list$pair_13
  pair_23 = result_list$pair_23
  
  pair_12.plot <- ggplot(pair_12, aes(y = rho, x = fct_rev(type), fill=type)) + 
    geom_boxplot() + ylab(TeX("$\\rho$")) +
    ylim(-1, 1) + geom_hline(yintercept = true_rho$rho12, size=2, color="red", linetype='dashed')+
    theme(axis.text.x = element_text(color = "grey20", size = 20),
          axis.text.y = element_text(color = "grey20", size = 20),  
          axis.title.x = element_text(color = "grey20", size = 0),
          axis.title.y = element_text(color = "grey20", size = 15),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15))
  
  pair_13.plot <- ggplot(pair_13, aes(y = rho, x = fct_rev(type), fill=type)) + 
    geom_boxplot() + ylab(TeX("$\\rho$")) +
    ylim(-1, 1) + geom_hline(yintercept = true_rho$rho13, size=2, color="red", linetype='dashed')+
    theme(axis.text.x = element_text(color = "grey20", size = 20),
          axis.text.y = element_text(color = "grey20", size = 20),  
          axis.title.x = element_text(color = "grey20", size = 0),
          axis.title.y = element_text(color = "grey20", size = 15),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15))
  
  pair_23.plot <- ggplot(pair_23, aes(y = rho, x = fct_rev(type), fill=type)) + 
    geom_boxplot() + ylab(TeX("$\\rho$")) + 
    ylim(-1, 1) + geom_hline(yintercept = true_rho$rho23, size=2, color="red", linetype='dashed')+
    theme(axis.text.x = element_text(color = "grey20", size = 20),
          axis.text.y = element_text(color = "grey20", size = 20),  
          axis.title.x = element_text(color = "grey20", size = 0),
          axis.title.y = element_text(color = "grey20", size = 15),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15))
  ggarrange(pair_12.plot, pair_13.plot, pair_23.plot, legend = "none", ncol=3)
} 

hist_1_pair <- function(pair, name, rho) {
  ggplot(pair, aes(x = rho, y = fct_rev(type), fill=type, height = stat(density)   )) + 
    geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE) +
    xlim(-1, 1)+geom_vline(xintercept = rho, size=0.75, color="red")+
    labs(title=name)+
    theme(axis.text.x = element_text(color = "grey20", size = 20),
          axis.text.y = element_text(color = "grey20", size = 20),  
          axis.title.x = element_text(color = "grey20", size = 15),
          axis.title.y = element_text(color = "grey20", size = 0),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15))
} 


summary_3_region <- function(summary_list, rho_vec, ...) {
  pair_12 = summary_list$pair_12
  pair_13 = summary_list$pair_13
  pair_23 = summary_list$pair_23
  
  pair_12_summary <- cbind(rho.true=rho_vec$rho12,
                           pair_12 %>% group_by(type) %>% summarise(RMSE=sqrt(mean((rho - rho_vec$rho12)^2)),
                                                                    BiasAbs = abs(mean(rho) - rho_vec$rho12),
                                                                    SD = sd(rho)))
  pair_13_summary <- cbind(rho.true=rho_vec$rho13,
                           pair_13 %>% group_by(type) %>% summarise(RMSE=sqrt(mean((rho - rho_vec$rho13)^2)),
                                                                    BiasAbs = abs(mean(rho) - rho_vec$rho13),
                                                                    SD = sd(rho)))
  pair_23_summary <- cbind(rho.true=rho_vec$rho23,
                           pair_23 %>% group_by(type) %>% summarise(RMSE=sqrt(mean((rho - rho_vec$rho23)^2)),
                                                                    BiasAbs = abs(mean(rho) - rho_vec$rho23),
                                                                    SD = sd(rho)))
  df <- rbind(pair_12_summary, pair_13_summary, pair_23_summary)
  
  #print(xtable(df,  digits = 4), include.rownames=FALSE, digits = 4)
}
