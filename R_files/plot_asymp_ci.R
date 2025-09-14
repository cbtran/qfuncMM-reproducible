library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
ggthemr::ggthemr("fresh")

coverage_results <- read_csv("out/std/asymp_ci_coverage_detailed.csv", col_types = "ffidfddfld", show_col_types = FALSE)
coverage_results$delta <- forcats::fct_relevel(coverage_results$delta, "low", "mid", "high")
coverage_results$psi <- forcats::fct_relevel(coverage_results$psi, "low", "mid", "high")

p <- coverage_results |>
  filter(method %in% c("reml")) |>
  group_by(delta, psi, level, rho_true) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = rho_true, linetype = rho_true, shape = rho_true
  )) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_line(position = position_dodge(width = 0.01)) +
  geom_errorbar(aes(
    ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
  ), width = 0.02, position = position_dodge(width = 0.01)) +
  geom_point(position = position_dodge(width = 0.01)) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  scale_color_brewer(palette = "Set2", labels = c("r12", "r13", "r23")) +
  scale_linetype_discrete(labels = c("r12", "r13", "r23")) +
  scale_shape_discrete(labels = c("r12", "r13", "r23")) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_detailed_reml.png", p, width = 10, height = 8)

p <- coverage_results |>
  filter(method %in% c("vecchia")) |>
  group_by(delta, psi, level, rho_true) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = rho_true, linetype = rho_true, shape = rho_true
  )) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_line(position = position_dodge(width = 0.01)) +
  geom_errorbar(aes(
    ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
  ), width = 0.02, position = position_dodge(width = 0.01)) +
  geom_point(position = position_dodge(width = 0.01)) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  scale_color_brewer(palette = "Set2", labels = c("r12", "r13", "r23")) +
  scale_linetype_discrete(labels = c("r12", "r13", "r23")) +
  scale_shape_discrete(labels = c("r12", "r13", "r23")) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_detailed_vecchia.png", p, width = 10, height = 8)

p <- coverage_results |>
  filter(method %in% c("reml", "ca", "ca_adjusted")) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, group = method
  )) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_line(position = position_dodge(width = 0.01)) +
  geom_errorbar(aes(
    ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
  ), width = 0.02, position = position_dodge(width = 0.01)) +
  geom_point(position = position_dodge(width = 0.01)) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  # theme_few() +
  scale_color_brewer(palette = "Set2", labels = c("ReML", "CA", "CA adjusted")) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_with_ca.png", p, width = 10, height = 8)

p <- coverage_results |>
  filter(method %in% c("vecchia", "vecchia_approx")) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, group = method
  )) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_line(position = position_dodge(width = 0.01)) +
  geom_errorbar(aes(
    ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
  ), width = 0.02, position = position_dodge(width = 0.01)) +
  geom_point(position = position_dodge(width = 0.01)) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  # theme_few() +
  scale_color_brewer(palette = "Set2") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_vecchia_approx.png", p, width = 10, height = 8)

p <- coverage_results |>
  filter(method %in% c("reml", "reml_approx")) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, group = method
  )) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_line(position = position_dodge(width = 0.01)) +
  geom_errorbar(aes(
    ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
  ), width = 0.02, position = position_dodge(width = 0.01)) +
  geom_point(position = position_dodge(width = 0.01)) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  # theme_few() +
  scale_color_brewer(palette = "Set2") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_approx.png", p, width = 10, height = 8)
