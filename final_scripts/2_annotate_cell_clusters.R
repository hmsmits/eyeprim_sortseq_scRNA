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

seurat_obj <- readRDS("/Users/hsmits7/surfdrive/CTI/eyeprim/final_data/merged_eyeprim_data_sct.rds")

# seurat_obj <- PrepSCTFindMarkers(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, assay = "SCT")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, algorithm = 4)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

png("./final_figures/total_umap.png")
DimPlot(seurat_obj)
dev.off()

clust_mark <- FindAllMarkers(seurat_obj, group.by = "seurat_clusters", features = VariableFeatures(seurat_obj[["SCT"]]))
top10 <- clust_mark %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

basal_markers <- c(
  "S100A2", "MIR205HG", "KRT14", "DST", "ROBO1", "COL17A1", "IGFBP2", "LGALS1", "SNHG29", "ATP1B3",
  "GAPDH", "SERPINF1", "TP63", "KRT17", "IGFBP7", "TP53I3P1", "GPNMB", "MTA2", "ZFP36L2",
  "NDUFA4L2", "CXCL14", "MT1E", "ARHGAP24", "IL18", "IFI1TM3", "IGFBP6", "KRT15", "KITLG"
)

keratinocyte_markers <- c(
  "AQP5", "LYPD2", "PSCA", "LCN2", "S100A8", "S100A9", "SLPI", "S100A4", "CXCL17", "CLU",
  "SERPINB4", "LGALS3", "FAM3D", "F3", "UPK3BL1", "SLC34A2", "MDK", "AGR2", "CLDN7", "WFDC2",
  "IFIT2", "FABP3", "SAT1", "CEACAM6", "FTH1", "TSPAN1", "RARRES1", "UPK1B",
  "CST3", "CD55", "BAG1", "TNFSF10", "MS4A1", "TMSB4X", "GSTA1", "TACSTD2", "C15orf48", "CYBA",
  "ELF3", "AC023154.2", "FAM3B", "SMM22", "VSIR", "KRT4", "KRT7", "AP001267.3", "VAMP5", "CLINT1",
  "CLIC3", "PLAAT4", "S100A11", "S100P", "B3GALT5-AS1", "ATP10B", "HAZF1", "ADGRF1", "CEACAM5",
  "CHFR", "PDZK1IP1", "PPPDEF", "CYP4B1", "RHOV", "STX3", "C9orf16", "ASAH1", "KRT13", "SPINK5",
  "A4GALT", "CTSS", "LINC01889", "BCAS1", "PERP", "B4GALT5", "GPX2", "GSTP1", "AQP3", "MUC20"
)

png("./final_figures/top_10_markers.png")
DoHeatmap(seurat_obj, top10$gene) + ggtitle("Top 10 markers")
dev.off()

png("./final_figures/basal_cell_markers.png")
DoHeatmap(seurat_obj, basal_markers) + ggtitle("Basal cell markers")
dev.off()

png("./final_figures/keratinocytes_cell_markers.png")
DoHeatmap(seurat_obj, keratinocyte_markers) + ggtitle("Keratinocyte cell markers")
dev.off()

# Based on the keratoncyte markers and basal cell markers we can say that cluster 1 are kerats
# And cluster 2 are basal cell
# Cluster 5 and 6 are immune cell/fibroblasts, in which we're not interested
# Cluster 3 and 4 are more difficult to distinguish 

seurat_obj <- subset(seurat_obj, seurat_clusters %in% c(1,2,3,4))
# Rename cluster IDs to biologically relevant cell type labels
new_cluster_names <- c(
  "1" = "keratinocytes", 
  "2" = "basal cells", 
  "3" = "clst_3", 
  "4" = "clst_4"
)

seurat_obj <- SetIdent(seurat_obj, value = "seurat_clusters")
seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)
seurat_obj$celltype <- Idents(seurat_obj)

saveRDS(seurat_obj, "./final_data/filtered_seurat_obj_sct.rds")
write.csv(top10, "./final_results/markers_seurat_clusters.csv")
