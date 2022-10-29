suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(ggridges))

output_dir <- "simulation-output/"
true_rho <- c(0.1, 0.35, 0.6)
args <- commandArgs(trailingOnly = TRUE)

M <- as.integer(args[1])
num_sims <- as.integer(args[2])

####### Boxplot #######
########## Read data function ########################
read_results <- function(result_folder) {
  pair.pearson <- read_csv(paste0(result_folder, "/rho_pearson.csv"))[-1]
  pair.inter <- read_csv(paste0(result_folder, "/rho_inter.csv"))[-1]
  
  pair.12 <- data.frame(type=factor(c(rep(c("REML", "FE", "CA"), each=num_sims)),
                                    levels=c("REML", "FE", "CA")),
                        rho=c(pair.inter$rho_12, pair.pearson$cor_intra_12, pair.pearson$cor_vanilla_12))
  
  pair.13 <- data.frame(type=factor(c(rep(c("REML", "FE", "CA"), each=num_sims)),
                                    levels=c("REML", "FE", "CA")),
                        rho=c(pair.inter$rho_13, pair.pearson$cor_intra_13, pair.pearson$cor_vanilla_13))
  
  pair.23 <- data.frame(type=factor(c(rep(c("REML", "FE", "CA"), each=num_sims)),
                                    levels=c("REML", "FE", "CA")),
                        rho=c(pair.inter$rho_23, pair.pearson$cor_intra_23, pair.pearson$cor_vanilla_23))
  
  pair.all <- cbind(rbind(pair.12, pair.13, pair.23),
                    inter=factor(c(rep(as.character(true_rho), each=3*num_sims)),
                                 levels=as.character(true_rho)))
  list(pair_12=pair.12, pair_13=pair.13, pair_23=pair.23, pair_all=pair.all)
}

########## Plot function ########################
create_boxplot <- function(result_list, true_rho) {
  
  pair_12 = cbind(result_list$pair_12, Pair="Pair 1-2")
  pair_13 = cbind(result_list$pair_13, Pair="Pair 1-3")
  pair_23 = cbind(result_list$pair_23, Pair="Pair 2-3")
  rbind(pair_12, pair_13, pair_23)
} 


################### Weak signal ###############3

# weak signal - weak intra
weak_signal.weak_intra <- read_results(paste0(output_dir, "weak-signal-weak-intra"))
weak_weak <- cbind(create_boxplot(weak_signal.weak_intra, "weak_signal_weak_intra"), signal ="weak", intra="weak")

# weak signal - strong intra
weak_signal.strong_intra <- read_results(paste0(output_dir, "weak-signal-strong-intra"))
weak_strong <- cbind(create_boxplot(weak_signal.strong_intra, "weak_signal_strong_intra"), signal = "weak", intra="strong")

################### Med signal ###############3

# med signal - weak intra
med_signal.weak_intra <- read_results(paste0(output_dir, "med-signal-weak-intra"))
med_weak <- cbind(create_boxplot(med_signal.weak_intra, "med_signal_weak_intra"), signal = "strong", intra="weak")

# med signal - strong intra
med_signal.strong_intra <- read_results(paste0(output_dir, "med-signal-strong-intra"))
med_strong <- cbind(create_boxplot(med_signal.strong_intra, "med_signal_strong_intra"), signal = "stron", intra="strong")

################### Strong signal ###############3

# strong signal - weak intra
strong_signal.weak_intra <- read_results(paste0(output_dir, "strong-signal-weak-intra"))
strong_weak <- cbind(create_boxplot(strong_signal.weak_intra, "strong_signal_weak_intra"), signal = "verystrong", intra="weak")

# strong signal - strong intra
strong_signal.strong_intra <- read_results(paste0(output_dir, "strong-signal-strong-intra"))
strong_strong <- cbind(create_boxplot(strong_signal.strong_intra, "strong_signal_strong_intra"), signal = "verystrong", intra="strong")



boxplot_df <- rbind(med_weak, med_strong, weak_weak, weak_strong, strong_weak, strong_strong)
boxplot_df$signal <- factor(boxplot_df$signal, levels=c("weak",  "strong", "verystrong"))
boxplot_df$intra <- factor(boxplot_df$intra, levels=c("weak",  "strong"))

levels(boxplot_df$signal) <- c(weak = TeX("$k_\\eta = 0.5$"),
                               strong = TeX("$k_\\eta = 1$"),
                               verystrong = TeX("$k_\\eta = 1.5$"))
levels(boxplot_df$intra) <- c(weak = TeX("$\\phi_{\\gamma} = 1$"),
                              strong = TeX("$\\phi_{\\gamma} = 0.25$"))

boxplot_pair12 <- boxplot_df[boxplot_df$Pair == "Pair 1-2",]
boxplot_pair13 <- boxplot_df[boxplot_df$Pair == "Pair 1-3",]
boxplot_pair23 <- boxplot_df[boxplot_df$Pair == "Pair 2-3",]

fig1 <- ggplot(boxplot_pair12, aes(y = rho, x = fct_rev(type), fill=fct_rev(type))) + 
  geom_boxplot() + ylab(TeX("$\\rho$")) +
  facet_grid(vars(signal), vars(intra), labeller=label_parsed) +
  ylim(-1, 1) + geom_hline(yintercept = true_rho[1], size=1, color="red", linetype='dashed')+
  theme_bw() +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x = element_text(color = "grey20", size = 30),
        axis.text.y = element_text(color = "grey20", size = 30),
        axis.title.x = element_text(color = "grey20", size = 0),
        axis.title.y = element_text(color = "grey20", size = 40),
        legend.text.align = 0,
        legend.title = element_text(size = 40),
        legend.text=element_text(size=15),
        aspect.ratio=1,
        legend.position = "none",
        legend.key.size = grid::unit(2, "lines"),
        strip.text = element_text(size = 40))
fig1
pdf(paste0(output_dir, "pair_12.pdf"), width = 12, height = 12)
fig1
dev.off()


fig2 <- ggplot(boxplot_pair13, aes(y = rho, x = fct_rev(type), fill=fct_rev(type))) + 
  geom_boxplot() + ylab(TeX("$\\rho$")) +
  facet_grid(vars(signal), vars(intra), labeller=label_parsed) +
  ylim(-1, 1) + geom_hline(yintercept = true_rho[2], size=1, color="red", linetype='dashed')+
  theme_bw() +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x = element_text(color = "grey20", size = 30),
        axis.text.y = element_text(color = "grey20", size = 30),
        axis.title.x = element_text(color = "grey20", size = 0),
        axis.title.y = element_text(color = "grey20", size = 40),
        legend.text.align = 0,
        legend.title = element_text(size = 40),
        legend.text=element_text(size=15),
        aspect.ratio=1,
        legend.position = "none",
        legend.key.size = grid::unit(2, "lines"),
        strip.text = element_text(size = 40))
fig2
pdf(paste0(output_dir, "pair_13.pdf"), width = 12, height = 12)
fig2
dev.off()


fig3 <- ggplot(boxplot_pair23, aes(y = rho, x = fct_rev(type), fill=fct_rev(type))) + 
  geom_boxplot() + ylab(TeX("$\\rho$")) +
  facet_grid(vars(signal), vars(intra), labeller=label_parsed) +
  ylim(-1, 1) + geom_hline(yintercept = true_rho[3], size=1, color="red", linetype='dashed')+
  theme_bw() +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x = element_text(color = "grey20", size = 30),
        axis.text.y = element_text(color = "grey20", size = 30),
        axis.title.x = element_text(color = "grey20", size = 0),
        axis.title.y = element_text(color = "grey20", size = 40),
        legend.text.align = 0,
        legend.title = element_text(size = 40),
        legend.text=element_text(size=15),
        aspect.ratio=1,
        legend.position = "none",
        legend.key.size = grid::unit(2, "lines"),
        strip.text = element_text(size = 40))
fig3
pdf(paste0(output_dir, "pair_23.pdf"), width = 12, height = 12)
fig3
dev.off()





####### Confidence interval #######
########## Read data function ########################
read_results <- function(result_folder, alpha) {
  
  # REML
  asymp_var_inter <- read.csv(paste0(result_folder, "/asymp_var_inter.csv"))[,-1]
  rho_inter <- read.csv(paste0(result_folder, "/rho_inter.csv"))[,-1]
  rho_inter_transformed <- atanh(rho_inter)
  
  reml_pair.12 <- data.frame(rho = rho_inter$rho_12,
                             lower = tanh(rho_inter_transformed$rho_12 - 1/sqrt(asymp_var_inter$var_12)*qnorm(1 - alpha/2)),
                             upper = tanh(rho_inter_transformed$rho_12 + 1/sqrt(asymp_var_inter$var_12)*qnorm(1 - alpha/2)))
  reml_pair.13 <- data.frame(rho = rho_inter$rho_13,
                             lower = tanh(rho_inter_transformed$rho_13 - 1/sqrt(asymp_var_inter$var_13)*qnorm(1 - alpha/2)),
                             upper = tanh(rho_inter_transformed$rho_13 + 1/sqrt(asymp_var_inter$var_13)*qnorm(1 - alpha/2)))
  reml_pair.23 <- data.frame(rho = rho_inter$rho_23,
                             lower = tanh(rho_inter_transformed$rho_23 - 1/sqrt(asymp_var_inter$var_23)*qnorm(1 - alpha/2)),
                             upper = tanh(rho_inter_transformed$rho_23 + 1/sqrt(asymp_var_inter$var_23)*qnorm(1 - alpha/2)))
  
  
  # Pearson
  pearson_inter <- read.csv(paste0(result_folder, "/rho_pearson.csv"))[,-1]
  pearson_inter_transformed <- atanh(pearson_inter)
  
  pearson_pair.12 <- data.frame(rho = pearson_inter$cor_vanilla_12,
                                lower = tanh(pearson_inter_transformed$cor_vanilla_12 - 1/sqrt(M-3)*qnorm(1 - alpha/2)),
                                upper = tanh(pearson_inter_transformed$cor_vanilla_12 + 1/sqrt(M-3)*qnorm(1 - alpha/2)))
  pearson_pair.13 <- data.frame(rho = pearson_inter$cor_vanilla_13,
                                lower = tanh(pearson_inter_transformed$cor_vanilla_13 - 1/sqrt(M-3)*qnorm(1 - alpha/2)),
                                upper = tanh(pearson_inter_transformed$cor_vanilla_13 + 1/sqrt(M-3)*qnorm(1 - alpha/2)))
  pearson_pair.23 <- data.frame(rho = pearson_inter$cor_vanilla_23,
                                lower = tanh(pearson_inter_transformed$cor_vanilla_23 - 1/sqrt(M-3)*qnorm(1 - alpha/2)),
                                upper = tanh(pearson_inter_transformed$cor_vanilla_23 + 1/sqrt(M-3)*qnorm(1 - alpha/2)))
  coverage <- data.frame(reml_12 = mean(reml_pair.12$upper >= true_rho[1] & reml_pair.12$lower <= true_rho[1]),
                         reml_13 = mean(reml_pair.13$upper >= true_rho[2] & reml_pair.13$lower <= true_rho[2]),
                         reml_23 = mean(reml_pair.23$upper >= true_rho[3] & reml_pair.23$lower <= true_rho[3]),
                         pearson_12 = mean(pearson_pair.12$upper >= true_rho[1] & pearson_pair.12$lower <= true_rho[1]),
                         pearson_13 = mean(pearson_pair.13$upper >= true_rho[2] & pearson_pair.13$lower <= true_rho[2]),
                         pearson_23 = mean(pearson_pair.23$upper >= true_rho[3] & pearson_pair.23$lower <= true_rho[3]))
  width <- data.frame(reml_12 = mean(abs(atanh(reml_pair.12$upper) - atanh(reml_pair.12$lower))),
                      reml_13 = mean(abs(atanh(reml_pair.13$upper) - atanh(reml_pair.13$lower))),
                      reml_23 = mean(abs(atanh(reml_pair.23$upper) - atanh(reml_pair.23$lower))),
                      pearson_12 = mean(abs(atanh(pearson_pair.12$upper) - atanh(pearson_pair.12$lower))),
                      pearson_13 = mean(abs(atanh(pearson_pair.13$upper) - atanh(pearson_pair.13$lower))),
                      pearson_23 = mean(abs(atanh(pearson_pair.23$upper) - atanh(pearson_pair.23$lower))))
  
  
  list(reml_pair_12=reml_pair.12, 
       reml_pair_13=reml_pair.13, 
       reml_pair_23=reml_pair.23,
       pearson_pair_12=pearson_pair.12,
       pearson_pair_13=pearson_pair.13,
       pearson_pair_23=pearson_pair.23,
       coverage=coverage,
       width=width)
}

create_CI <- function(result_folder) {
  med_med_coverage <- c()
  med_med_width <- c()
  for (alpha in c(0.05, 0.1, 0.15, 0.2, 0.25)) {
    test = med_signal.med_intra <- read_results(result_folder, alpha)
    med_med_coverage <- rbind(med_med_coverage, cbind(alpha=1-alpha, test$coverage))
    med_med_width <-  rbind(med_med_width, cbind(alpha=1-alpha, test$width))
  }
  
  coverage_df <- data.frame(alpha=rep(med_med_coverage$alpha, 3),
                            coverage = c(med_med_coverage$reml_12,
                                         med_med_coverage$reml_13,
                                         med_med_coverage$reml_23),
                            Pair = rep(c("Pair 1-2", "Pair 1-3", "Pair 2-3"), each=5))
  coverage_df
  
}


# Create df

weak_weak = cbind(create_CI(paste0(output_dir, "weak-signal-weak-intra")), signal ="weak", intra="weak")
weak_strong =  cbind(create_CI(paste0(output_dir, "weak-signal-strong-intra")), signal = "weak", intra="strong")
med_weak =  cbind(create_CI(paste0(output_dir, "med-signal-weak-intra")), signal = "strong", intra="weak")
med_strong =  cbind(create_CI(paste0(output_dir, "med-signal-strong-intra")), signal = "strong", intra="strong")
strong_weak =  cbind(create_CI(paste0(output_dir, "strong-signal-weak-intra")), signal = "verystrong", intra="weak")
strong_strong =  cbind(create_CI(paste0(output_dir, "strong-signal-strong-intra")), signal = "verystrong", intra="strong")

CI_df <- rbind(med_weak, med_strong, weak_weak, weak_strong, strong_weak, strong_strong)
write.csv(CI_df, 'CI_all.csv')
CI_df$signal <- factor(CI_df$signal, levels=c("weak",  "strong", "verystrong"))
CI_df$intra <- factor(CI_df$intra, levels=c("weak",  "strong"))

levels(CI_df$signal) <- c(weak = TeX("$k_\\eta = 0.5$"),
                          strong = TeX("$k_\\eta = 1$"),
                          verystrong = TeX("$k_\\eta = 1.5$"))
levels(CI_df$intra) <- c(weak = TeX("$\\phi_{\\gamma} = 1$"),
                         strong = TeX("$\\phi_{\\gamma} = 0.25$"))

# Plot
fig <-   ggplot(CI_df, aes(x=alpha, y=coverage, color=Pair, linetype=Pair, shape=Pair)) +
  geom_point(size=2.5) +
  geom_line(size=1) +
  labs(y="Coverage", x=TeX("$1 - \\alpha$")) +
  scale_linetype_manual(values = c("dotdash", "dashed", "longdash")) +
  facet_grid(vars(signal), vars(intra), labeller=label_parsed) +
  # facet_wrap(~type, scales = "free_x", 
  #            labeller=label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(color = "grey20", size = 15),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        legend.text.align = 0,
        legend.title = element_text(size=20),
        legend.text=element_text(size=15),
        aspect.ratio=1,
        legend.position = "right",
        legend.key.size = grid::unit(2, "lines"),
        strip.text = element_text(size=20)) +
  geom_abline(intercept = 0, slope = 1)
fig
pdf(paste0(output_dir,"reml_ci.pdf"), width = 12, height = 12)
fig
dev.off()

# Summary table

weak_weak = cbind(k_eta=0.5, phi_gamma=1, read.csv(paste0(output_dir, "weak-signal-weak-intra/weak-signal-weak-intra.csv")))
weak_strong =  cbind(k_eta=0.5, phi_gamma=0.25, read.csv(paste0(output_dir, "weak-signal-strong-intra/weak-signal-strong-intra.csv")))
med_weak =  cbind(k_eta=1, phi_gamma=1, read.csv(paste0(output_dir, "med-signal-weak-intra/med-signal-weak-intra.csv")))
med_strong =  cbind(k_eta=1, phi_gamma=0.25, read.csv(paste0(output_dir, "med-signal-strong-intra/med-signal-strong-intra.csv")))
strong_weak =  cbind(k_eta=1.5, phi_gamma=1, read.csv(paste0(output_dir, "strong-signal-weak-intra/strong-signal-weak-intra.csv")))
strong_strong =  cbind(k_eta=1.5, phi_gamma=0.25, read.csv(paste0(output_dir, "strong-signal-strong-intra/strong-signal-strong-intra.csv")))
write.csv(rbind(weak_weak, weak_strong, med_weak, med_strong, strong_weak, strong_strong), 'results_all.csv')


