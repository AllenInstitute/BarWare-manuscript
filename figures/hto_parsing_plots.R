library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)

res_df <- fread("BarCounter/Pool-32-HTO_Tag_Counts.csv")


p1 <- ggplot() +
  geom_histogram(
    data = res_df,
    aes(x = HT2),
    bins = 100,
    fill = "darkgreen"
  ) +
  scale_x_continuous(
    "HTO UMIs",
    limits = c(0, 2000),
    expand = c(0,0)) +
  scale_y_continuous(
    "N Barcodes",
    limits = c(0, 1500),
    expand = c(0,0)
  ) +
  theme_bw(base_size = 6)
  
p2 <- ggplot() +
  geom_histogram(
    data = res_df %>%
      filter(HT2 > 10),
    aes(x = HT2),
    bins = 100,
    fill = "darkgreen"
  ) +
  scale_x_continuous(
    "HTO UMIs",
    limits = c(0, 2000),
    expand = c(0,0)) +
  scale_y_continuous(
    "N Barcodes",
    limits = c(0, 1500),
    expand = c(0,0)
  ) +
  theme_bw(base_size = 6)
  
p3 <- ggplot() +
  geom_histogram(
    data = res_df %>%
      filter(HT2 > 10),
    aes(x = log10(HT2)),
    bins = 90,
    fill = "darkgreen"
  ) +
  scale_x_continuous(
    "log10(HTO UMIs)",
    limits = c(1, 3.6),
    expand = c(0,0)) +
  scale_y_continuous(
    "N Barcodes",
    expand = c(0,0)
  ) +
  theme_bw(base_size = 6)
  
vals <- log10(res_df$HT2[res_df$HT2 > 10])

km <- kmeans(vals, 2)

km_centers <- data.frame(
  center_pos = km$centers[,1]
)

center_delta <- data.frame(
  log_diff = max(km_centers$center_pos) - min(km_centers$center_pos),
  fold_diff = 10^max(km_centers$center_pos) / 10^min(km_centers$center_pos),
  xpos = sum(km_centers$center_pos) / 2
)

km_df <- data.frame(
  val = 10^vals,
  cl = km$cluster
)

p4 <- ggplot() +
  geom_histogram(
    data = km_df,
    aes(x = log10(val),
        fill = as.factor(cl)),
    bins = 90
  ) +
  geom_vline(
    data = km_centers,
    aes(xintercept = center_pos),
    color = "black"
  ) +
  geom_text(
    data = center_delta,
    aes(x = xpos,
        y = 1000,
        label = round(fold_diff, 2))
  ) +
  scale_x_continuous(
    "log10(HTO UMIs)",
    limits = c(1, 3.6),
    expand = c(0,0)) +
  scale_y_continuous(
    "N Barcodes",
    expand = c(0,0)
  ) +
  theme_bw(base_size = 6) +
  theme(legend.position = "none")

high_clust <- which(km$centers == max(km$centers))

cutoff <- data.frame(
  xpos = min(vals[km$cluster == high_clust]),
  val = 10^min(vals[km$cluster == high_clust])
)

p5 <- ggplot() +
  geom_histogram(
    data = km_df,
    aes(x = log10(val),
        fill = as.factor(cl)),
    bins = 90
  ) +
  geom_vline(
    data = cutoff,
    aes(xintercept = xpos),
    color = "dodgerblue"
  ) +
  geom_text(
    data = cutoff,
    aes(x = xpos,
        y = 1000,
        label = round(val,2)),
    hjust = -0.1
  ) +
  scale_x_continuous(
    "log10(HTO UMIs)",
    limits = c(1, 3.6),
    expand = c(0,0)) +
  scale_y_continuous(
    "N Barcodes",
    expand = c(0,0)
  ) +
  theme_bw(base_size = 6) +
  theme(legend.position = "none")

joint_plot <- plot_grid(
  p1, p2, p3, p4, p5,
  nrow = 1
)

ggsave(
  "HTO_parsing_histograms.pdf",
  joint_plot,
  width = 7.5,
  height = 0.7
)
