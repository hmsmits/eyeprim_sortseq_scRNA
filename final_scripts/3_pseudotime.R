library(monocle3)
library(Seurat)
library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library("flowCore")
library(readxl)
library(edgeR)
library(tidyverse)
library(SeuratWrappers)

filter <- dplyr::filter
select <- dplyr::select

source("/Users/hsmits7/surfdrive/CTI/eyeprim/scripts/plot_3d_cells.R")

path <- "/Users/hsmits7/surfdrive/CTI/eyeprim/"
setwd(path)

seurat_obj <- readRDS("/Users/hsmits7/surfdrive/CTI/eyeprim/final_data/filtered_seurat_obj_sct.rds")
seurat_obj <- RunPCA(seurat_obj, assay = "SCT")
ElbowPlot(seurat_obj, ndims = 50)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50)

cds <- as.cell_data_set(seurat_obj)
cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = c("leiden"))

cluster <- seurat_obj$celltype
names(cluster) <- cds@colData@rownames
cluster <- as.factor(cluster)
cds@clusters$UMAP$clusters <- cluster 	

cds <- learn_graph(cds)
cds <- order_cells(cds)

png("./final_figures/pseudo_umap_clusters_seurat.png",  height = 2000, width = 3000, res = 300)
DimPlot(seurat_obj, group.by = "celltype")
dev.off()

png("./final_figures/pseudo_umap_clusters.png",  height = 2000, width = 3000, res = 300)
plot_cells(cds, color_cells_by = "cluster")
dev.off()

png("./final_figures/pseudotime_umap_.png",  height = 2000, width = 3000, res = 300)
plot_cells(cds,  color_cells_by = "pseudotime", cell_size = 0.9,scale_to_range = T, )
dev.off()

# graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=4)
# 
# saveRDS(graph_test_res, "./final_results/pseudo_time_results.rds")
graph_test_res <- readRDS("./final_results/pseudo_time_results.rds")

plotgenes <- graph_test_res %>%  filter(q_value < 0.05) %>% arrange(q_value, 1/abs(morans_test_statistic), decreasing = T) %>% rownames()

antimicrobial_genes <- c(
  "PIGR", "WFDC2", "S100A4", "S100A8", 
  "S100A9", "SLPI", "LCN2", "CXCL17",
  "LTF",  "PIGR"
)

AQPs <- plotgenes[grepl("AQP", plotgenes)]
S100 <- plotgenes[grepl("S100", plotgenes)]

png("./final_figures/top20_pseudogenes.png", height = 3000, width = 3000, res = 300)
FeaturePlot(seurat_obj, plotgenes[1:20])
dev.off()

png("./final_figures/anti_microbial_genes.png", height = 3000, width = 3000, res = 300)
FeaturePlot(seurat_obj,antimicrobial_genes)
dev.off()

png("./final_figures/AQP_genes.png", height = 1200, width = 3000,   res = 300)
FeaturePlot(seurat_obj,AQPs)
dev.off()

png("./final_figures/S100_genes.png", height = 2000, width = 3000, res = 300)
FeaturePlot(seurat_obj, S100)
dev.off()

