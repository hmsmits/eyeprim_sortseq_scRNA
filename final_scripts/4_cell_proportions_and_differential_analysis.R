# library(monocle3)
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
# library(SeuratWrappers)
library(clusterProfiler)
library(org.Hs.eg.db)  # Annotation database for human genes
library(enrichplot)  # For visualization
library(openxlsx)

perform_GSEA <- function(celltype_list){
  GSEA_res <- list()
  for(celltype_name in names(celltype_list)){
    celltype <- celltype_list[[celltype_name]]
    celltype_sort <- celltype %>% arrange(-logFC)
    gene_list <- celltype_sort$logFC
    names(gene_list) <- celltype_sort$gene
    
    Entrez_IDs <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop = F)  
    message("number of duplicated symbols: ", sum(duplicated(Entrez_IDs$SYMBOL)))
    Entrez_IDs <- Entrez_IDs %>% filter(!duplicated(SYMBOL))
    Entrez_IDs <- Entrez_IDs[match(Entrez_IDs$SYMBOL, names(gene_list)),]
    Entrez_list <- celltype_sort$logFC
    names(Entrez_list) <- Entrez_IDs$ENTREZID
    Entrez_list <- Entrez_list[!is.na(names(Entrez_list))]
    
    gsea_kegg <- gseKEGG(geneList = Entrez_list, 
                         organism = "hsa", 
                         # nPerm = 1000, 
                         # minGSSize = 10, 
                         maxGSSize = 500, 
                         pvalueCutoff = 0.05)
    
    gsea_go <- gseGO(geneList = Entrez_list, 
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                     # nPerm = 1000, 
                     # minGSSize = 10, 
                     maxGSSize = 500, 
                     pvalueCutoff = 0.05)
    
    GSEA_res_celltype <- list(KEGG = gsea_kegg, GO = gsea_go)
    GSEA_res[[celltype_name]] <- GSEA_res_celltype
  }
  return(GSEA_res)
}

# Function to perform differential expression using edgeR's LRT
perform_EdgeR_LRT <- function(srt_obj, cluster, DE_list, contrast, paired = F){
  # Subset object for specific cell type
  sub <- subset(srt_obj, celltype == cluster)
  # Get the aggregated expression for the pseudo bulk, assay is the raw RNA assay because we will perform EdgeR default normalization!!!
  pseudo <- AggregateExpression(sub, group.by = c("timepoint", "SampleID", "PatientID"), return.seurat = T, assays = "RNA")
  
  if(paired){
    design = model.matrix(~ timepoint + PatientID, data = pseudo@meta.data)
  }else{
    design = model.matrix(~ timepoint, data = pseudo@meta.data)
  }
  
  # edgeR pipeline
  y = DGEList(counts = pseudo@assays$RNA$counts, group =  pseudo$timepoint)
  y = calcNormFactors(y, method = "TMM")
  y = estimateGLMRobustDisp(y, design)
  fit = glmFit(y, design = design)
  test = glmLRT(fit, coef = 2)
  
  res <- topTags(test, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')  %>% 
    mutate(cell_type = cluster, coef = colnames(design)[2])
  
  DE_list[[cluster]] <- res
  return(DE_list)
}

filter <- dplyr::filter
select <- dplyr::select

source("/Users/hsmits7/surfdrive/CTI/eyeprim/scripts/plot_3d_cells.R")

path <- "/Users/hsmits7/surfdrive/CTI/eyeprim/"
setwd(path)

seurat_obj <- readRDS("/Users/hsmits7/surfdrive/CTI/eyeprim/eyeprim_sortseq_scRNA/final_data/filtered_seurat_obj_sct.rds")

cluster_counts <- table(
  Cluster = seurat_obj$celltype,
  # Sample = combined_seurat_sct$SampleID,
  Treatment = seurat_obj$Treament,
  Timepoint = seurat_obj$timepoint
)

cluster_df <- as.data.frame(cluster_counts) %>% 
  mutate(Cluster = factor(Cluster, levels = c("keratinocytes", "basal cells", "clst_3", "clst_4")))

# Calculate percentages
cluster_df <- cluster_df %>%
  group_by(Timepoint, Treatment) %>%
  mutate(Percent = Freq / sum(Freq) * 100)
cluster_df <- cluster_df %>% filter(!Freq == 0)
cluster_df <- cluster_df %>% mutate(Timepoint = factor(Timepoint, levels = c("Pso", "BL", "4weeks")), Treatment = factor(Treatment, levels = c("Pso", "Dupi", "Tralo")))
# Plot the data

total_percentages <- ggplot(cluster_df, aes(x = Timepoint, y = Percent, fill = Cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Treatment, scales = "free") +
  ylab("Percentage of Cells") +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggsci::scale_fill_bmj()

total_absolute <- ggplot(cluster_df, aes(x = Timepoint, y = Freq, fill = Cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Treatment, scales = "free_x") +
  ylab("Total number of Cells") +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggsci::scale_fill_bmj()

celltype_absolute <- ggplot(cluster_df, aes(x = Timepoint, y = Freq, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Cluster, scales = "free_x") +
  ylab("Total number of Cells") +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggsci::scale_fill_d3()

png("./eyeprim_sortseq_scRNA/final_figures/cellproportions_total.png")
total_percentages
dev.off()

png("./eyeprim_sortseq_scRNA/final_figures/cell_absolute_total.png")
total_absolute
dev.off()

png("./eyeprim_sortseq_scRNA/final_figures/celltype_absolute_total.png")
celltype_absolute
dev.off()

cluster_counts <- table(
  Cluster = seurat_obj$celltype,
  Patient = seurat_obj$PatientID,
  Treatment = seurat_obj$Treament,
  Timepoint = seurat_obj$timepoint,
  Sample = seurat_obj$SampleID
)

cluster_df <- as.data.frame(cluster_counts) %>% 
  mutate(Cluster = factor(Cluster, levels = c("keratinocytes", "basal cells", "clst_3", "clst_4")))

# Calculate percentages
cluster_df <- cluster_df %>%
  group_by(Timepoint, Treatment, Sample) %>%
  mutate(Percent = Freq / sum(Freq) * 100)
cluster_df <- cluster_df %>% filter(!Freq == 0)
cluster_df <- cluster_df %>% mutate(Timepoint = factor(Timepoint, levels = c("Pso", "BL", "4weeks")), Treatment = factor(Treatment, levels = c("Pso", "Dupi", "Tralo")))
# Plot the data

# List of treatments to loop through
treatments <- c("Pso", "Dupi", "Tralo")

for (treatment in treatments) {
  
  # Percentage Plot
  p_percentage <- cluster_df %>%
    filter(Treatment == treatment) %>%
    ggplot(aes(x = Timepoint, y = Percent, fill = Cluster)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Treatment + Patient, scales = "free") +
    ylab("Percentage of Cells") +
    xlab("Treatment") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggsci::scale_fill_bmj() +
    scale_x_discrete(labels = c("BL" = "Baseline", "4W" = "4 Weeks"))
  
  # Save Percentage Plot
  ggsave(filename = paste0("./eyeprim_sortseq_scRNA/final_figures/", treatment, "_percentage.png"),
         plot = p_percentage, 
         width = 10, height = 6, dpi = 300)
  
  # Absolute Plot
  p_absolute <- cluster_df %>%
    filter(Treatment == treatment) %>%
    ggplot(aes(x = Timepoint, y = Freq, fill = Cluster)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Treatment + Patient, scales = "free_x") +
    ylab("Absolute Cell Count") +
    xlab("Treatment") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggsci::scale_fill_bmj() +
    scale_x_discrete(labels = c("BL" = "Baseline", "4W" = "4 Weeks"))
  
  # Save Absolute Plot
  ggsave(filename = paste0("./eyeprim_sortseq_scRNA/final_figures/", treatment, "_absolute.png"),
         plot = p_absolute, 
         width = 10, height = 6, dpi = 300)
}

DE_markers <- list()

for(treatment in treatments){
  if(treatment == "Pso"){
    low_counts_patients <- c("control 1") # has low counts
    seurat_subs <- subset(seurat_obj, timepoint %in%c("Pso", "BL") & !PatientID%in% low_counts_patients)
    paired <- F
  }else{
    low_counts_patients <- c("dupi 3", "tralo 10") # both follow up patients have low counts
    seurat_subs <- subset(seurat_obj, Treament == treatment & !PatientID%in% low_counts_patients)
    paired <- T
  }
  patients <- unique(seurat_subs$PatientID)
  message("running pseudo Edger LRT with the following samples: ", paste0(patients, collapse = ", "))
  DE_markers_sublist <- list()
  for(cluster in unique(seurat_subs$celltype)){
    try(DE_markers_sublist <- perform_EdgeR_LRT(seurat_subs, cluster = cluster, DE_list = DE_markers_sublist, paired = paired))
  }
  DE_markers[[treatment]] <- DE_markers_sublist
}

for(treatment in treatments){
  DE_markers_trt <- DE_markers[[treatment]]
  if(treatment == "Pso"){
    title <- "AD vs Pso"
    file <- "vlc_plot_AD_vs_Pso.png"
  }else{
    title <- paste0(treatment, " BL vs. 4 weeks of treatment")
    file <- paste0("vlc_plot_", treatment, "_BL_vs_4weeks.png")
    
  }
  DE_markers_trt_pseudo_df <- Reduce(bind_rows, DE_markers_trt) %>% mutate(label = case_when(FDR < 0.05 ~ gene, T ~ NA))
  volcano <- DE_markers_trt_pseudo_df %>% mutate(Gene = case_when(FDR < 0.05 ~ gene, T ~ NA)) %>% 
    ggplot(aes(x = logFC, y = -log(FDR), col = FDR < 0.05, label = Gene)) + 
    geom_point() + theme_bw() + ggrepel::geom_label_repel(max.overlaps = 5) + 
    ggtitle(title) + facet_wrap(~cell_type, scales = "free")
  
  ggsave(filename = paste0("./eyeprim_sortseq_scRNA/final_figures/", file),
         plot = volcano, 
         width = 10, height = 6, dpi = 300)
}

saveRDS(DE_markers, "./eyeprim_sortseq_scRNA/final_results/DE_markers.rds")

for(treatment in treatments){
  DE_markers_trt <- DE_markers[[treatment]] 
  wb <- createWorkbook()
  
  file <- paste0("./eyeprim_sortseq_scRNA/final_results/", treatment, "_DE_genes_EdgeR_LRT.xlsx")
  
  # Add each data frame to a separate sheet
  for (celltype in names(DE_markers_trt)) {
    addWorksheet(wb, celltype)
    writeData(wb, celltype, DE_markers_trt[[celltype]])
  }
  
  # Save the workbook
  saveWorkbook(wb, file, overwrite = TRUE)
}


GSEA_list <- list()

for(treatment in treatments){
  DE_markers_trt <- DE_markers[[treatment]]
  if(treatment == "Pso"){
    title_GO <- "GO: AD vs Pso"
    title_KEGG <- "KEGG: AD vs Pso"
    file_GO <- "GO_plot_AD_vs_Pso.png"
    file_KEGG <- "KEGG_plot_AD_vs_Pso.png"
  }else{
    title_GO <- paste0("GO: ", treatment, " BL vs. 4 weeks of treatment")
    title_KEGG <- paste0("KEGG: ", treatment, " BL vs. 4 weeks of treatment")
    file_GO <- paste0("GO_plot_", treatment, "_BL_vs_4weeks.png")
    file_KEGG <- paste0("KEGG_plot_", treatment, "_BL_vs_4weeks.png")
    
  }
  GSEA_list[[treatment]] <- perform_GSEA(DE_markers_trt)
  
  for(celltype in names(GSEA_list[[treatment]])){
    GSEA_trt <- GSEA_list[[treatment]]
    GSEA_celltype <- GSEA_trt[[celltype]]
    if(length(GSEA_celltype$GO$Description) > 0){
      dtplt_GO <- dotplot(GSEA_celltype$GO) + ggtitle(paste0(celltype, ": ", title_GO))
      ggsave(filename = paste0("./eyeprim_sortseq_scRNA/final_figures/", paste0(celltype, "_",file_GO)),
             plot = dtplt_GO, 
             width = 6, height = 10, dpi = 300)
    }
    if(length(GSEA_celltype$KEGG$Description) > 0){
      dtplt_KEGG <- dotplot(GSEA_celltype$KEGG) + ggtitle(paste0(celltype, ": ", title_KEGG))
      ggsave(filename = paste0("./eyeprim_sortseq_scRNA/final_figures/", paste0(celltype, "_",file_KEGG)),
             plot = dtplt_KEGG, 
             width = 6, height = 10, dpi = 300)
      
    }
  }
}


# seurat_obj$patient_timepoint <- paste(seurat_obj$PatientID, seurat_obj$timepoint, sep = "_")
# 
# library(patchwork)
# 
# plots <- FeaturePlot(
#   subset(seurat_obj, celltype == "basal cells" & Treament == "Tralo"),
#   features = "KRT14",
#   split.by = "patient_timepoint", cols = 2, combine = F
# ) + patchwork::plot_layout(ncol = 2)
# 
# RidgePlot(subset(seurat_obj, celltype == "basal cells" & Treament == "Tralo" & KRT14 > 0), features = "KRT14", ncol = 2)
# RidgePlot(subset(seurat_obj, celltype %in% c("basal cells" & KRT14 > 0), features = "KRT14", ncol = 2)
# 
# # Set identities to that
# Idents(seurat_obj) <- "patient_timepoint"
# 
# # Plot
# VlnPlot(seurat_obj, features = "KRT14")
# FeaturePlot(subset(seurat_obj, celltype == "basal cells" & Treament == "Tralo"), "KRT14", split.by = "timepoint")
