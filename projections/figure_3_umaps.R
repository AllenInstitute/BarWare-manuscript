library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(BarMixer)
library(cowplot)
library(immutils)
library(ArchR)
library(pheatmap)

setwd("/mnt/disks/barcode-tender-manuscript/cell-hashing-pipeline/analysis")

# read in H5s
h5_dir <- "/mnt/disks/barcode-tender-manuscript/cell-hashing-pipeline/merge_h5/"
mem1 <- read_h5_seurat(file.path(h5_dir, "X017-P1_2735BW-MEM-1.h5"))
mem2 <- read_h5_seurat(file.path(h5_dir, "X017-P1_2735BW-MEM-2.h5"))
niv1 <- read_h5_seurat(file.path(h5_dir, "X017-P1_2735BW-NIV-1.h5"))
niv2 <- read_h5_seurat(file.path(h5_dir, "X017-P1_2735BW-NIV-2.h5"))
non1 <- read_h5_seurat(file.path(h5_dir, "X017-P1_2735BW-NON-1.h5"))
non2 <- read_h5_seurat(file.path(h5_dir, "X017-P1_2735BW-NON-2.h5"))

# merge Seurat objects
merged <- merge(mem1, y = c(mem2, niv1, niv2, non1, non2), project = "bartender")
merged


# RNA QC -------------------------------
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        group.by = 'pbmc_sample_id')
plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged)
plot1 + plot2

# evaluate mito percentage
merged@meta.data %>% group_by(pbmc_sample_id) %>%
  summarise(min_mito = min(percent.mt),
            max_mito = max(percent.mt),
            median_mito = median(percent.mt),
            mean_mito = mean(percent.mt),
            iqr_mito = IQR(percent.mt),
            cutoff = (sum(percent.mt > 25) / length(percent.mt))*100)

# evaluate umi counts
merged@meta.data %>% group_by(pbmc_sample_id) %>%
  summarise(min_umi = min(nCount_RNA),
            median_umi = median(nCount_RNA),
            mean_umi = mean(nCount_RNA),
            max_umi = max(nCount_RNA))
merged@meta.data %>% group_by(pbmc_sample_id) %>%
  summarise(perc = (sum(nCount_RNA > 25000) / length(nCount_RNA))*100)
merged@meta.data %>% group_by(pbmc_sample_id) %>%
  summarise(perc = (sum(nCount_RNA < 1000) / length(nCount_RNA))*100)

# evaluate genes detected
merged@meta.data %>% group_by(pbmc_sample_id) %>%
  summarise(min_feature = min(nFeature_RNA),
            median_feature = median(nFeature_RNA),
            mean_feature = mean(nFeature_RNA),
            max_feature = max(nFeature_RNA))
merged@meta.data %>% group_by(pbmc_sample_id) %>%
  summarise(perc = (sum(nFeature_RNA < 500) / length(nFeature_RNA))*100)


# filter cells with high mito, low genes detected, or outlier UMI counts ----------------------------------

# originally 113,414 singlets
merged <- subset(merged, subset = percent.mt < 25 & nCount_RNA < 25000 & nCount_RNA >= 1000 & nFeature_RNA >= 500)
# 105,395 singlets afterwards
merged

# add metadata field for FACS population
condition <- sub("-[1|2]","",merged@meta.data$pbmc_sample_id)
condition <- sub("2735BW-","", condition)
condition[condition == "MEM"] <- "Memory T"
condition[condition == "NIV"] <- "Naive T"
condition[condition == "NON"] <- "Non-T"

merged@meta.data$population <- condition

# label transfer -----------------

# use Seurat PBMC CITE-seq reference
reference <- LoadH5Seurat("/shared/Seurat_References/pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2",
        label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# normalize x017 data with SCTransform to match reference
merged <- SCTransform(merged)
merged <- RunPCA(merged, assay = "SCT", reduction.name = 'sct_pca')
merged <- RunUMAP(merged, assay = 'SCT', reduction = 'sct_pca', dims = 1:50)
merged <- FindNeighbors(object = merged, dims = 1:50, reduction = 'sct_pca')
merged <- FindClusters(object = merged)

# find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  reference.assay = "SCT",
  query = merged,
  recompute.residuals = FALSE,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# transfer cell type labels and protein data from the reference to the query
merged <- MapQuery(
  anchorset = anchors,
  query = merged,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)


# create plots --------------------

# create a cluster confusion matrix across each population using the confusionMatrix() function

# subset l2 to remove labels with < 50 cells
table(merged$predicted.celltype.l2)
low_l2 <- c("ASDC","CD4 Proliferating","pDC","cDC1","Eryth","CD8 Proliferating","Platelet")
subset_l2 <- subset(x = merged, predicted.celltype.l2 %ni% low_l2) # removed 137 cells

cM <- confusionMatrix(paste0(subset_l2$predicted.celltype.l2), paste0(subset_l2$population))
# plot confusion matrix as a heatmap
cM <- cM / Matrix::rowSums(cM)
clusterHeatmap <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  cellwidth = 70,
  cellheight = 20
)


# get hexadecimal values for viridis palette
library(scales)
show_col(viridis_pal(option = "C")(20))

umap_sample <- DimPlot(merged, reduction = "umap",
                       group.by = "pbmc_sample_id",
                       raster = FALSE) +
  labs(title = "PBMC Replicates")
umap_sample

umap_pop <- DimPlot(merged,reduction = "umap",
                    group.by = "population",
                    raster = FALSE,
                    cols = c("Memory T" = "#20068FFF","Naive T" = "#7AD7F0","Non-T" = "#D6556DFF")) +
  guides(color = guide_legend(override.aes = list(size=12), ncol=1))
umap_pop

l1_umap = DimPlot(merged, reduction = "umap",
                  group.by = "predicted.celltype.l1",
                  raster = FALSE,
                  label = FALSE,
                  repel = TRUE, legend.title.size = 10) +
  labs(title = "Seurat Celltype Labels (level 1)")
l1_umap


# marker genes
markerGenes  <- c("CD3D", "MS4A1","CD14","GNLY") # MS4A1 = CD20, FCGR3A = CD16
naive_mem_markers <- c("CD3D","IL7R","CCR7","S100A4")


# red = #D2042D
# purple = #0D0887FF

marker_umap <- FeaturePlot(subset_l2,
                           features = markerGenes,
                           min.cutoff = "q10",
                           raster = FALSE,
                           cols = c("grey","#D2042D"))
marker_umap

naive_mem_umap <- FeaturePlot(subset_l2, features = naive_mem_markers,
                              min.cutoff = "q10", cols = c("grey","#0D0887FF"),
                              raster = FALSE)
naive_mem_umap



# IGKC, IGLC1 and IGLC2
# B cell light chain genes
light_chain <- c("IGKC","IGLC1","IGLC2")
FeaturePlot(subset_l2, features = light_chain, min.cutoff = "q10")


# save Seurat object ---------------------
saveRDS(merged, file = "fig3_merged_x017_so.rds")
merged <- readRDS(file = "fig3_merged_x017_so.rds")

