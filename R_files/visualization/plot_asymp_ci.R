library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(patchwork)
ggthemr::ggthemr("fresh")

coverage_results <- read_csv("out/std/asymp_ci_coverage_detailed.csv", col_types = "ffidfddfld", show_col_types = FALSE)
coverage_results$delta <- forcats::fct_relevel(coverage_results$delta, "low", "mid", "high")
coverage_results$psi <- forcats::fct_relevel(coverage_results$psi, "low", "mid", "high")

# Define consistent visual scheme for methods
method_colors <- c(
  "reml" = "#1f77b4",
  "reml_approx" = "#aec7e8",
  "vecchia" = "#ff7f0e",
  "vecchia_approx" = "#ffbb78",
  "ca" = "#2ca02c", # "#72b764",
  "ca_adjusted" = "#2ca02c"
)

method_linetypes <- c(
  "reml" = "solid",
  "reml_approx" = "dotted",
  "vecchia" = "dashed",
  "vecchia_approx" = "dotted",
  "ca" = "dashed", # "dotted",
  "ca_adjusted" = "dashed"
)

method_shapes <- c(
  "reml" = 16,
  "reml_approx" = 1,
  "vecchia" = 17,
  "vecchia_approx" = 2,
  "ca" = 0,
  "ca_adjusted" = 15
)

method_labels <- c(
  "reml" = "ReML",
  "reml_approx" = "ReML Approx",
  "vecchia" = "Vecchia",
  "vecchia_approx" = "Vecchia Approx",
  "ca" = "CA",
  "ca_adjusted" = "CA Adjusted"
)

# Define method groups to avoid repetition
method_groups <- list(
  "reml_ca" = c("reml", "ca", "ca_adjusted"),
  "vecchia_approx" = c("vecchia", "vecchia_approx"),
  "reml_approx" = c("reml", "reml_approx"),
  "reml_vecchia" = c("reml", "vecchia")
)

# Function to get scales for a method group
get_method_scales <- function(methods) {
  list(
    color = scale_color_manual(values = method_colors[methods], labels = method_labels[methods]),
    linetype = scale_linetype_manual(values = method_linetypes[methods], labels = method_labels[methods]),
    shape = scale_shape_manual(values = method_shapes[methods], labels = method_labels[methods])
  )
}

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
ggsave("plots/asymp_ci_coverage_detailed_reml.pdf", p, width = 10, height = 8, device = cairo_pdf)

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
ggsave("plots/asymp_ci_coverage_detailed_vecchia.pdf", p, width = 10, height = 8, device = cairo_pdf)

p <- coverage_results |>
  filter(method %in% method_groups$reml_ca) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, linetype = method, shape = method, group = method
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
  get_method_scales(method_groups$reml_ca)$color +
  get_method_scales(method_groups$reml_ca)$linetype +
  get_method_scales(method_groups$reml_ca)$shape +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_with_ca.pdf", p, width = 10, height = 8, device = cairo_pdf)

p <- coverage_results |>
  filter(method %in% method_groups$vecchia_approx) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, linetype = method, shape = method, group = method
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
  get_method_scales(method_groups$vecchia_approx)$color +
  get_method_scales(method_groups$vecchia_approx)$linetype +
  get_method_scales(method_groups$vecchia_approx)$shape +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_vecchia_approx.pdf", p, width = 10, height = 8, device = cairo_pdf)

p <- coverage_results |>
  filter(method %in% method_groups$reml_approx) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, linetype = method, shape = method, group = method
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
  get_method_scales(method_groups$reml_approx)$color +
  get_method_scales(method_groups$reml_approx)$linetype +
  get_method_scales(method_groups$reml_approx)$shape +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_approx.pdf", p, width = 10, height = 8, device = cairo_pdf)

# Additional plot: psi=mid only, 1x2 facet for delta=low,high
p <- coverage_results |>
  filter(method %in% method_groups$reml_ca, psi == "mid", delta %in% c("low", "high")) |>
  group_by(delta, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, linetype = method, shape = method, group = method
  )) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_line(position = position_dodge(width = 0.01)) +
  geom_errorbar(aes(
    ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
  ), width = 0.02, position = position_dodge(width = 0.01)) +
  geom_point(position = position_dodge(width = 0.01)) +
  facet_wrap(~delta,
    ncol = 2,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x)
    )
  ) +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  get_method_scales(method_groups$reml_ca)$color +
  get_method_scales(method_groups$reml_ca)$linetype +
  get_method_scales(method_groups$reml_ca)$shape +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5)
  )
p
ggsave("plots/asymp_ci_coverage_psi_mid.pdf", p, width = 8, height = 4, device = cairo_pdf)

# 9x9 grid plot: ReML vs Vecchia across all delta-psi combinations
p <- coverage_results |>
  filter(method %in% method_groups$reml_vecchia) |>
  group_by(delta, psi, level, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, linetype = method, shape = method, group = method
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
  get_method_scales(method_groups$reml_vecchia)$color +
  get_method_scales(method_groups$reml_vecchia)$linetype +
  get_method_scales(method_groups$reml_vecchia)$shape +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )
p
ggsave("plots/asymp_ci_coverage_reml_vs_vecchia.pdf", p, width = 12, height = 10, device = cairo_pdf)

# 2x2 plot using patchwork: specific delta-psi combinations (without CA Adjusted)
# Function to create individual plot
create_plot <- function(delta_val, psi_val) {
  coverage_results |>
    filter(method %in% c("reml", "ca"), delta == delta_val, psi == psi_val) |>
    group_by(level, method) |>
    summarize(
      coverage = mean(covered, na.rm = TRUE),
      coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
    ) |>
    ggplot(aes(
      x = as.numeric(as.character(level)), y = coverage,
      color = method, linetype = method, shape = method, group = method
    )) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_line(position = position_dodge(width = 0.01)) +
    geom_errorbar(aes(
      ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
    ), width = 0.02, position = position_dodge(width = 0.01)) +
    geom_point(position = position_dodge(width = 0.01)) +
    labs(
      x = "1 - \u03B1",
      y = "Coverage",
      title = paste0("\u03B4: ", delta_val, ", \u03C8: ", psi_val)
    ) +
    scale_color_manual(values = method_colors[c("reml", "ca")], labels = method_labels[c("reml", "ca")]) +
    scale_linetype_manual(values = method_linetypes[c("reml", "ca")], labels = method_labels[c("reml", "ca")]) +
    scale_shape_manual(values = method_shapes[c("reml", "ca")], labels = method_labels[c("reml", "ca")]) +
    scale_y_continuous(limits = c(0.2, 1)) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 10)
    )
}

# Create individual plots
p1 <- create_plot("low", "mid")
p2 <- create_plot("high", "mid")
p3 <- create_plot("mid", "low")
p4 <- create_plot("mid", "high")

# Combine plots using patchwork
combined_plot <- (p1 | p2) / (p3 | p4)

# Add a shared legend at the bottom
combined_plot <- combined_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

combined_plot
ggsave("plots/asymp_ci_coverage_2x2_facet.pdf", combined_plot, width = 10, height = 8, device = cairo_pdf)

# 2x2 plot using patchwork: ReML, CA, and Vecchia methods
# Function to create individual plot with 3 methods
create_plot_3methods <- function(delta_val, psi_val) {
  coverage_results |>
    filter(method %in% c("reml", "ca", "vecchia"), delta == delta_val, psi == psi_val) |>
    group_by(level, method) |>
    summarize(
      coverage = mean(covered, na.rm = TRUE),
      coverage_sd = 1.96 * sd(covered, na.rm = TRUE) / sqrt(n())
    ) |>
    ggplot(aes(
      x = as.numeric(as.character(level)), y = coverage,
      color = method, linetype = method, shape = method, group = method
    )) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_line(position = position_dodge(width = 0.01)) +
    geom_errorbar(aes(
      ymin = coverage - coverage_sd, ymax = coverage + coverage_sd
    ), width = 0.02, position = position_dodge(width = 0.01)) +
    geom_point(position = position_dodge(width = 0.01)) +
    labs(
      x = "1 - \u03B1",
      y = "Coverage",
      title = paste0("\u03B4: ", delta_val, ", \u03C8: ", psi_val)
    ) +
    scale_color_manual(values = method_colors[c("reml", "ca", "vecchia")], labels = method_labels[c("reml", "ca", "vecchia")]) +
    scale_linetype_manual(values = method_linetypes[c("reml", "ca", "vecchia")], labels = method_labels[c("reml", "ca", "vecchia")]) +
    scale_shape_manual(values = method_shapes[c("reml", "ca", "vecchia")], labels = method_labels[c("reml", "ca", "vecchia")]) +
    scale_y_continuous(limits = c(0.2, 1)) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 10)
    )
}

# Create individual plots with 3 methods
p1_3m <- create_plot_3methods("low", "mid")
p2_3m <- create_plot_3methods("high", "mid")
p3_3m <- create_plot_3methods("mid", "low")
p4_3m <- create_plot_3methods("mid", "high")

# Combine plots using patchwork
combined_plot_3m <- (p1_3m | p2_3m) / (p3_3m | p4_3m)

# Add a shared legend at the bottom
combined_plot_3m <- combined_plot_3m +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

combined_plot_3m
ggsave("plots/asymp_ci_coverage_2x2_facet_3methods.pdf", combined_plot_3m, width = 10, height = 8, device = cairo_pdf)
