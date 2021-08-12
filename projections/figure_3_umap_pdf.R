library(ggplot2)
library(ggrastr)
library(pheatmap)
library(dplyr)
library(immutils)
library(BarMixer)

umap_df <- read.csv("figure_3_umap_data.csv.gz")
umap_df$replicate <- sub("2735BW-","",umap_df$pbmc_sample_id)

l2_table <- table(umap_df$predicted.celltype.l2)
low_l2 <- names(l2_table)[l2_table <= 100]
umap_df$predicted.celltype.l2[umap_df$predicted.celltype.l2 %in% low_l2] <- "Other"

umap_df <- umap_df %>%
  mutate(cell_class = case_when(
    predicted.celltype.l2 %in% c("CD4 Naive","CD8 Naive","dnT") ~ "Naive T",
    grepl("^CD[4|8]",predicted.celltype.l2) ~ "Memory T",
    predicted.celltype.l2 %in% c("MAIT","Treg","gdT") ~ "Memory T",
    TRUE ~ "Non-T"
  ))

class_summary <- umap_df %>%
  group_by(population) %>%
  mutate(n_cells = n()) %>%
  ungroup() %>%
  group_by(population, cell_class) %>%
  summarise(n_class = n(),
            frac_population = n() / n_cells[1],
            perc_population = round(100 * n() / n_cells[1], 2))
  

cM <- table(umap_df$predicted.celltype.l2, umap_df$population)
cM <- as.matrix(cM)

# plot confusion matrix as a heatmap
cM <- cM / rowSums(cM)

archr_whiteblue <- c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")

clusterHeatmap <- pheatmap(
  mat = cM, 
  color = colorRampPalette(archr_whiteblue)(256),
  border_color = "black",
  cellwidth = 70,
  cellheight = 20
)

pdf("confusion_matrix.pdf", width = 6, height = 10)
clusterHeatmap
dev.off()

umap_df <- umap_df[sample(1:nrow(umap_df), nrow(umap_df), replace = F),]

# FACS gating plots

facs_umap <- ggplot() +
  geom_point_rast(
    data = umap_df,
    aes(x = umap_1,
        y = umap_2,
        color = population),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  scale_color_manual("FACS Population",
                     breaks = c("Memory T","Naive T","Non-T"),
                     values = c("#20068F","#7AD7F0","#D6556D")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.height = unit(1, "inches"))

ggsave("facs_umap.pdf",
       facs_umap,
       width = 5, height = 6)

facs_umap_nolegend <- ggplot() +
  geom_point_rast(
    data = umap_df,
    aes(x = umap_1,
        y = umap_2,
        color = population),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  scale_color_manual("FACS Population",
                     breaks = c("Memory T","Naive T","Non-T"),
                     values = c("#20068F","#7AD7F0","#D6556D")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none")

ggsave("facs_umap_nolegend.pdf",
       facs_umap_nolegend,
       width = 5, height = 5)

# Replicate plots

rep_umap <- ggplot() +
  geom_point_rast(
    data = umap_df,
    aes(x = umap_1,
        y = umap_2,
        color = replicate),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  scale_color_manual("Replicate",
                     breaks = c("MEM-1","MEM-2",
                                "NIV-1","NIV-2",
                                "NON-1","NON-2"),
                     values = c("#20068F","#796ABC",
                                "#7AD7F0","#3D6C78",
                                "#D6556D","#E699A7")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.5, "inches"))

ggsave("rep_umap.pdf",
       rep_umap,
       width = 5, height = 6)

rep_umap_nolegend <- ggplot() +
  geom_point_rast(
    data = umap_df,
    aes(x = umap_1,
        y = umap_2,
        color = replicate),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  scale_color_manual("Replicate",
                     breaks = c("MEM-1","MEM-2",
                                "NIV-1","NIV-2",
                                "NON-1","NON-2"),
                     values = c("#20068F","#796ABC",
                                "#7AD7F0","#3D6C78",
                                "#D6556D","#E699A7")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none")

ggsave("rep_umap_nolegend.pdf",
       rep_umap_nolegend,
       width = 5, height = 5)

# Well plots

well_ids <- paste0("X017-P1C1W",3:8)
well_colors <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
well_umap <- ggplot() +
  geom_point_rast(
    data = umap_df[sample(1:nrow(umap_df), nrow(umap_df), replace = F),],
    aes(x = umap_1,
        y = umap_2,
        color = well_id),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  scale_color_manual("Well ID",
                     breaks = well_ids,
                     values = well_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.5, "inches"))

ggsave("well_umap.pdf",
       well_umap,
       width = 5, height = 6)

well_umap_nolegend <- ggplot() +
  geom_point_rast(
    data = umap_df[sample(1:nrow(umap_df), nrow(umap_df), replace = F),],
    aes(x = umap_1,
        y = umap_2,
        color = well_id),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  scale_color_manual("Well ID",
                     breaks = well_ids,
                     values = well_colors) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none")

ggsave("well_umap_nolegend.pdf",
       well_umap_nolegend,
       width = 5, height = 5)

# Cell Type umap

type_umap <- ggplot() +
  geom_point_rast(
    data = umap_df,
    aes(x = umap_1,
        y = umap_2,
        color = predicted.celltype.l2),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  theme_bw(base_size = 7) +
  scale_color_varibow() +
  large_guides() +
  theme(legend.position = "bottom")

ggsave("type_umap.pdf",
       type_umap,
       width = 5, height = 5)

type_umap_nolegend <- ggplot() +
  geom_point_rast(
    data = umap_df,
    aes(x = umap_1,
        y = umap_2,
        color = predicted.celltype.l2),
    raster.dpi = 300,
    size = 0.001,
    shape = 20
  ) +
  theme_bw(base_size = 7) +
  scale_color_varibow() +
  theme(legend.position = "none")

ggsave("type_umap_nolegend.pdf",
       type_umap_nolegend,
       width = 5, height = 5)

# well vs replicate barplot
well_distribution_barplot <- qc_aligned_barplot(
  umap_df,
  category_x = "replicate",
  name_x = "Replicate",
  category_y = "well_id",
  category_name = "Well ID",
  colorset_y = "rainbow"
)

ggsave("well_distribution_barplot.pdf",
       well_distribution_barplot,
       width = 5, height = 5)


## Marker plots

marker_df <- read.csv("figure3_marker_counts.csv.gz", row.names = 1)

marker_cols <- names(marker_df)[-1]

umap_df <- left_join(umap_df, marker_df)

norm_umap_df <- umap_df
for(marker_col in marker_cols) {
  norm_umap_df[[marker_col]] <- log10(norm_umap_df[[marker_col]] / norm_umap_df$n_umis * 1e4 + 1)
}

marker_plot_list <- lapply(
  marker_cols,
  function(marker) {
    pmarker <- rlang::parse_expr(marker)
    
    ggplot() +
      geom_point_rast(
        data = norm_umap_df %>% arrange(!!pmarker),
        aes(x = umap_1,
            y = umap_2,
            color = !!pmarker),
        raster.dpi = 300,
        size = 0.001,
        shape = 20
      ) +
      theme_bw(base_size = 7) +
      scale_color_gradient(low = "gray90", high = "red") +
      theme(legend.position = "bottom") +
      ggtitle(marker)
  }
)

lapply(
  1:length(marker_cols),
  function(i) {
    marker <- marker_cols[[i]]
    plot <- marker_plot_list[[i]]
    
    ggsave(paste0("marker_umap_",marker,".pdf"),
           plot,
           width = 5, height = 5.5)
  }
)

## Marker Dotplot

dotplot_markers <- rev(c("CD3D", "ITGB1","S100A4","CCR7","HLA.DRA","CST3","CD79A","KLRD1"))

marker_dotplot_df <- norm_umap_df %>%
  select(replicate, one_of(dotplot_markers)) %>%
  reshape2::melt(
    id.vars = c("replicate"),
    variable.name = "gene",
    value.name = "norm_expr"
  ) %>%
  group_by(replicate, gene) %>%
  summarize(mean_expr = mean(norm_expr),
            med = median(norm_expr),
            med_nz = median(norm_expr[norm_expr > 0]),
            frac = sum(norm_expr > 0) / n())

marker_dotplot <- ggplot() +
  geom_point(data = marker_dotplot_df,
             aes(x = replicate,
                 y = gene,
                 size = frac,
                 fill = med_nz),
             color = "black",
             pch = 21) +
  scale_fill_gradient(low = "gray90", high = "red") +
  scale_size_area(breaks = c(0.25,0.5,0.75,1),
                  limits = c(0,1)) +
  theme_bw(base_size = 7)

ggsave("marker_dotplot.pdf",
       marker_dotplot,
       width = 2.5,
       height = 2.5)
