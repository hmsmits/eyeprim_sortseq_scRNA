library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library("flowCore")
library(readxl)
library(edgeR)

# Function to append metadata to a Seurat object
append_meta_data <- function(seurat_obj,meta){
  seurat_obj@meta.data$SampleID <- gsub("-raw", "",gsub(".*-(AE)-s(\\d+)(.*)", "\\1\\2\\3", seurat_obj@meta.data$orig.ident))
  new_metadata <- dplyr::left_join(seurat_obj@meta.data, meta)
  rownames(new_metadata) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj,new_metadata)
  return(seurat_obj)
}


# Define working directory and data paths
path <- "/Users/hsmits7/surfdrive/CTI/eyeprim/"
setwd(path)

# Define directories
datadir_batch1 <- "/Users/hsmits7/surfdrive/CTI/eyeprim/data/batch_1/raw_count_tables/non_poisson_corrected/"
datadir_batch2 <- "/Users/hsmits7/surfdrive/CTI/eyeprim/data/batch_2/raw_count_tables/non_poisson_corrected/"
datadir_batch3 <- "/Users/hsmits7/surfdrive/CTI/eyeprim/data/batch_3/raw_count_tables/non_poisson_corrected/"
datadir_batch4 <- "/Users/hsmits7/surfdrive/CTI/eyeprim/data/batch_4/raw_count_tables/non_poisson_corrected/"
# Load metadata
metadata <- read.csv("/Users/hsmits7/surfdrive/CTI/eyeprim/data/metadata.csv", sep = ";") 
features <- read.csv(
  "/Users/hsmits7/surfdrive/CTI/eyeprim/data/batch_1/raw_count_tables/poisson_corrected/UMC-AE-s004-raw/features.tsv",
  sep = "\t", col.names = c("ENSG", "Gene_Symbol", "Type")
)

# Collect folder names for each batch
folders <- c(
  list.files(datadir_batch1, full.names = TRUE),
  list.files(datadir_batch2, full.names = TRUE),
  list.files(datadir_batch3, full.names = TRUE),
  list.files(datadir_batch4, full.names = TRUE)
)
folders <- folders[!grepl(".pdf", folders)]  # Remove non-data folders

# Initialize lists and read barcode layout, cell counts, FACS mappings
seurat_list <- list()
batch_index <- data.frame(batch = basename(folders), index = 1:length(folders))
barcodes <- read.table("/Users/hsmits7/surfdrive/CTI/eyeprim/data/cellseq_barcode_position.tsv", skip=1, col.names = c('Well position', 'Name', 'Barcode'))
cell_counts <- read.csv("/Users/hsmits7/surfdrive/CTI/eyeprim/data/number_of_cells.csv", header = F, sep = ";")
cell_counts$V1 <- gsub("filtered", "raw", cell_counts$V1)

# Load FACS mapping and FCS files
mapping <- read_excel("/Users/hsmits7/surfdrive/CTI/eyeprim/fcs_sc_mapping.xlsx", col_names = c("scID", "SampleID", "Treament", "timepoint", "FACS_file"))
mapping$scID <- gsub("filtered", "raw", mapping$scID)
mapping <- na.omit(mapping)
fcs_list <- lapply(mapping$FACS_file, function(file){read.FCS(file.path("/Volumes/LAB/lti/_Clinical/Group-van-Wijk/Hidde/AD_epigenetics/Successful sorts/", file))})
names(fcs_list) <- mapping$scID
# fcs_list <- list
merged_fcs <- list()

# Loop over folders to construct Seurat objects (special case for AE026 handled separately)
# For each folder, read the matrix, filter spike-ins and mt genes, add FACS data if available
# Create one Seurat object per plate, and store in list
for (folder in folders) {
  if(basename(folder) == "UMC-AE-s026-raw"){
    message("fill in this part, split samples")
    seurat_obj <- ReadMtx(
      mtx = file.path(folder, "matrix.mtx"),
      features = file.path(folder, "features.tsv"),
      cells = file.path(folder, "barcodes.tsv"),
      feature.column = 2
    )
    FACS_data_tralo_3 <- fcs_list[["UMC-AE-s026-raw_tralo_3"]]
    FACS_data_tralo_4 <- fcs_list[["UMC-AE-s026-raw_tralo_4"]]
    
    # Reorder the cells so they are in the order of A1 -> A24 -> B24 -> B1
    # Also filter out the spike-ins wells, and the empty wells
    ordering <- order(match(colnames(seurat_obj),barcodes$Barcode))
    # ordering_tralo_3 <- 
    nCells_tralo_3 <- nrow(exprs(FACS_data_tralo_3))
    nCells_tralo_4 <- nrow(exprs(FACS_data_tralo_4))
    
    tralo_3_cells <- c(1:168, 192:(192-(nCells_tralo_3 - 169)))
    tralo_4_cells <- setdiff(1:374, tralo_3_cells)[1:nCells_tralo_4]
    
    seurat_obj <- seurat_obj[,ordering]
    
    seurat_obj_tralo_3 <- seurat_obj[,tralo_3_cells]
    seurat_obj_tralo_4 <- seurat_obj[,tralo_4_cells]
    
    seurat_obj_tralo_3 <- CreateSeuratObject(counts = seurat_obj_tralo_3, project = paste0(basename(folder), "_tralo_3"))
    seurat_obj_tralo_4 <- CreateSeuratObject(counts = seurat_obj_tralo_4, project = paste0(basename(folder), "_tralo_4"))
    
    seurat_obj_tralo_3$plate <- folder
    seurat_obj_tralo_4$plate <- folder
    
    seurat_obj_tralo_3[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_tralo_3, pattern = "^MT-")
    seurat_obj_tralo_4[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_tralo_4, pattern = "^MT-")
    
    seurat_obj_tralo_3[["percent.spike"]] <- PercentageFeatureSet(seurat_obj_tralo_3, pattern = "^ERCC-")
    seurat_obj_tralo_4[["percent.spike"]] <- PercentageFeatureSet(seurat_obj_tralo_4, pattern = "^ERCC-")
    
    seurat_obj_tralo_3 <- subset(seurat_obj_tralo_3, features = setdiff(rownames(seurat_obj_tralo_3), grep("^ERCC-", rownames(seurat_obj_tralo_3), value = TRUE)))
    seurat_obj_tralo_4 <- subset(seurat_obj_tralo_4, features = setdiff(rownames(seurat_obj_tralo_4), grep("^ERCC-", rownames(seurat_obj_tralo_4), value = TRUE)))
    
    seurat_obj_tralo_3 <- subset(seurat_obj_tralo_3, features = setdiff(rownames(seurat_obj_tralo_3), grep("^MT-", rownames(seurat_obj_tralo_3), value = TRUE)))
    seurat_obj_tralo_4 <- subset(seurat_obj_tralo_4, features = setdiff(rownames(seurat_obj_tralo_4), grep("^MT-", rownames(seurat_obj_tralo_4), value = TRUE)))
    
    message("tralo_3: ", basename(folder))
    message("number of cells seurat = ", ncol(seurat_obj_tralo_3), "\n", "number of cells FACS = ", nCells_tralo_3)
    message("tralo_4: ", basename(folder))
    message("number of cells seurat = ", ncol(seurat_obj_tralo_4), "\n", "number of cells FACS = ", nCells_tralo_4)
    
    seurat_obj_tralo_3$CD45 <- exprs(FACS_data_tralo_3)[,CD45_channel]
    seurat_obj_tralo_4$CD45 <- exprs(FACS_data_tralo_4)[,CD45_channel]
    
    seurat_obj_tralo_3$EPCAM <- exprs(FACS_data_tralo_3)[,epcam_channel]
    seurat_obj_tralo_4$EPCAM <- exprs(FACS_data_tralo_4)[,epcam_channel]
    
    seurat_obj_tralo_3$CD45_level <- ifelse(seurat_obj_tralo_3$CD45 > 10**2.3, "CD45+", "CD45-")
    seurat_obj_tralo_4$CD45_level <- ifelse(seurat_obj_tralo_4$CD45 > 10**2.3, "CD45+", "CD45-")
    
    seurat_obj_tralo_3$EPCAM_level <- ifelse(seurat_obj_tralo_3$EPCAM > 10**3.1, "EPCAM+", "EPCAM-")
    seurat_obj_tralo_4$EPCAM_level <- ifelse(seurat_obj_tralo_4$EPCAM > 10**3.1, "EPCAM+", "EPCAM-")
    
    # seurat_obj_tralo_3 <- subset(seurat_obj_tralo_3, subset = percent.spike < 15 & nCount_RNA > 1000 & nFeature_RNA < 5000 )
    # seurat_obj_tralo_4 <- subset(seurat_obj_tralo_4, subset = percent.spike < 15 & nCount_RNA > 1000 & nFeature_RNA < 5000 )
    
    seurat_list[[paste0(basename(folder), "_tralo_3")]] <- seurat_obj_tralo_3
    seurat_list[[paste0(basename(folder), "_tralo_4")]] <- seurat_obj_tralo_4
    
  }else{
    # to do: identify how many cells to use, based on FACS or on the patient data file 
    seurat_obj <- ReadMtx(
      mtx = file.path(folder, "matrix.mtx"),
      features = file.path(folder, "features.tsv"),
      cells = file.path(folder, "barcodes.tsv"),
      feature.column = 2
    )
    FACS_data <- fcs_list[[basename(folder)]]
    
    # Reorder the cells so they are in the order of A1 -> A24 -> B24 -> B1
    # Also filter out the spike-ins wells, and the empty wells
    ordering <- order(match(colnames(seurat_obj),barcodes$Barcode))
    # nCells <- nrow(exprs(FACS_data))
    if(is.null(FACS_data)){
      message("skipping")
      nCells <- cell_counts[cell_counts$V1 == basename(folder),"V5"]
      if(!basename(folder) %in% cell_counts$V1){
        nCells <- 376 
      }
      ordering <- ordering[1:nCells]
      seurat_obj <- seurat_obj[,ordering]
      seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = basename(folder))
      
      seurat_obj$plate <- folder
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      seurat_obj[["percent.spike"]] <- PercentageFeatureSet(seurat_obj, pattern = "^ERCC-")
      seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), grep("^ERCC-", rownames(seurat_obj), value = TRUE)))
      seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), grep("^MT-", rownames(seurat_obj), value = TRUE)))
      seurat_obj$CD45 <- rep(NA, nCells)
      seurat_obj$EPCAM <-  rep(NA, nCells)
    }else{
      nCells <- nrow(exprs(FACS_data))
      ordering <- ordering[1:nCells]
      seurat_obj <- seurat_obj[,ordering]
      seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = basename(folder))
      
      seurat_obj$plate <- folder
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      seurat_obj[["percent.spike"]] <- PercentageFeatureSet(seurat_obj, pattern = "^ERCC-")
      seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), grep("^ERCC-", rownames(seurat_obj), value = TRUE)))
      seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), grep("^MT-", rownames(seurat_obj), value = TRUE)))
      
      message(basename(folder))
      message("number of cells seurat = ", ncol(seurat_obj), "\n", "number of cells FACS = ", nCells)
      
      ## Append the CD45 and EPCAM data (unnormalized), last two columns are always first CD45 then EPCAM
      # FACS_data <- merged_fcs[[basename(folder)]]
      # time_column <- which(colnames(FACS_data) == "Time")
      # FACS_data <- FACS_data[,-time_column]
      ncol(FACS_data)
      epcam_channel <- colnames(FACS_data)[which(pData(parameters(FACS_data))$desc%in%c("EPCAM PE","EpCAM-PE"))]
      CD45_channel <- colnames(FACS_data)[which(pData(parameters(FACS_data))$desc%in%c("CD45 AF700","CD45-AF700"))]
      
      seurat_obj$CD45 <- exprs(FACS_data)[,CD45_channel]
      seurat_obj$EPCAM <- exprs(FACS_data)[,epcam_channel]
      
      seurat_obj$CD45_level <- ifelse(seurat_obj$CD45 > 10**2.3, "CD45+", "CD45-")
      seurat_obj$EPCAM_level <- ifelse(seurat_obj$EPCAM > 10**3.1, "EPCAM+", "EPCAM-")
      
    }
    # Filter low quality cells
    # seurat_obj <- subset(seurat_obj, subset = percent.spike < 15 & nCount_RNA > 1000 & nFeature_RNA < 5000 )
    seurat_list[[basename(folder)]] <- seurat_obj
  }
}



# seurat_list_sct <- lapply(seurat_list, function(seurat_obj) {
#   SCTransform(
#     seurat_obj,
#     method = "glmGamPoi",
#     # variable.features.n = 10000,
#     return.only.var.genes = FALSE,
#     # vars.to.regress = c("percent.mt"),
#     verbose = FALSE,
#     min_cells = 0,
#     vst.flavor="v2"
#   )
# })

# Merge individual objects into one combined Seurat object
combined_seurat_merged <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)])

png("./final_figures/vln_plots.png")
QC_plot <- VlnPlot(combined_seurat_sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.spike"), ncol = 3)
dev.off()

png("./final_figures/scttr_plots.png")
FeatureScatter(combined_seurat_sct, "nCount_RNA", "nFeature_RNA")
dev.off()

combined_seurat_sct <- subset(combined_seurat_sct, subset = percent.spike < 15 & nCount_RNA > 1000 & nFeature_RNA < 5000 )

png("./final_figures/vln_plots_flt.png")
QC_plot <- VlnPlot(combined_seurat_sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.spike"), ncol = 3)
dev.off()

png("./final_figures/scttr_plots_flt.png")
FeatureScatter(combined_seurat_sct, "nCount_RNA", "nFeature_RNA")
dev.off()

combined_seurat_sct <-   SCTransform(
  combined_seurat_merged,
  method = "glmGamPoi",
  # variable.features.n = 10000,
  return.only.var.genes = FALSE,
  # vars.to.regress = c("percent.mt"),
  verbose = FALSE,
  min_cells = 0,
  vst.flavor="v2"
)


combined_seurat_sct <- append_meta_data(combined_seurat_sct, metadata)
combined_seurat_sct <- FindVariableFeatures(combined_seurat_sct)
# vrbl_features_sct <- SelectIntegrationFeatures(object.list = seurat_list_sct, nfeatures = 2000)
# VariableFeatures(combined_seurat_sct[["SCT"]]) <- vrbl_features_sct
combined_seurat_sct <- PrepSCTFindMarkers(combined_seurat_sct)

saveRDS(seurat_list, "./final_data/seurat_list.rds")
saveRDS(combined_seurat_sct, "./final_data/merged_eyeprim_data_sct.rds")