# Generates the supplemental plot showing robustness of Vecchia neighbor selection to the number of neighbors m.

library(ggplot2)
ggthemr::ggthemr("fresh")

neighbors <- read.csv("out/std/stage2_vecchia_neighbors/results_mid-mid.csv")
vecchia <- read.csv("out/std/stage2_vecchia/results_mid-mid.csv")
vecchia$m_seq <- 100

df <- rbind(
  neighbors[, c("simid", "m_seq", "region1_uniqid", "region2_uniqid", "rho")],
  vecchia[, c("simid", "m_seq", "region1_uniqid", "region2_uniqid", "rho")]
)

df$region_pair <- factor(paste0("r", df$region1_uniqid, df$region2_uniqid))
df$m_seq <- factor(df$m_seq)
df$yintercept <- ifelse(df$region_pair == "r12", 0,
  ifelse(df$region_pair == "r13", 0.35, 0.6)
)

p <- ggplot(df, aes(x = region_pair, y = rho, fill = m_seq)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_segment(
    aes(
      x = as.numeric(region_pair) - 0.5,
      xend = as.numeric(region_pair) + 0.5,
      y = yintercept,
      yend = yintercept
    ),
    lty = 2, color = "black"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Region pair", y = expression(rho), fill = "# neighbors") +
  theme(legend.position = "bottom")

cairo_pdf(file.path("plots", "vecchia_neighbors.pdf"), width = 8, height = 4, onefile = TRUE)
print(p)
dev.off()
