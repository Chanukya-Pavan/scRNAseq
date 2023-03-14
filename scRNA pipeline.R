
# Packages
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyverse)
library(RCurl)
library(ggplot2)
library(readr)
library(gridExtra)
library(harmony)

# List of seurat objects in list format and list object format
objects <- list(UCinc_UFM_L169_Lung_IPFDonor_14kDirectNuclei_cDNA, UCinc_UFM_L169_Lung_IPFDonor_35kSortedNuclei_cDNA,UCinc_UMich_IPF0165_Lung_IPFDonor_14kDirectNuclei_cDNA,UCinc_UMich_IPF0165_Lung_IPFDonor_35kSortedNuclei_cDNA,UCinc_UNC_KKCFFT044N_Lung_HealthyDonor_14kDirectNuclei_cDNA, UCinc_UNC_KKCFFT044N_Lung_HealthyDonor_35kSortedNuclei_cDNA,UCinc_UFM_L0171_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UFM_L0172_Lung_IPFDonor_FACSCounted15kNuclei_cDNA,UCinc_UNC_DD0025N_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_UNC_DD026P_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA,UCinc_UC_001_030_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UFM_IPF0175_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UFM_L0173_Lung_IPFDonor_NoSort15kNuclei_cDNA,UCinc_UFM_L201_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UMich_L182_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UNC_DD016N_Lung_HealthyDonor_NoSort15kNuclei_cDNA,UCinc_UNC_DD039N_Lung_HealthyDonor_NoSort15kNuclei_cDNA, UCinc_UNC_DD057N_Lung_HealthyDonor_NoSort15kNuclei_cDNA,UCinc_UC_001_007_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UC_001_021_Lung_IPFDonor_FACSCounted15kNuclei_cDNA,UCinc_UC_001_026_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UC_UC_1_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA,UCinc_UFM_L195_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UNC_DD003O_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA,UCinc_UNC_DD020P_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_CCHMC_HL_CCHMC_01_Lung_HealthyDonor_FACSCounted16kNuclei_cDNA,UCinc_UFM_L210_Lung_IPFDonor_FACSCounted16kNuclei_cDNA, UCinc_UMich_UMILD146_LLL0212_Lung_IPFDonor_FACSCounted16kNuclei_cDNA,UCinc_UMILD0145_L0211_Lung_IPFDonor_FACSCounted16kNuclei_cDNA)
flag <- list("UCinc_UFM_L169_Lung_IPFDonor_14kDirectNuclei_cDNA","UCinc_UFM_L169_Lung_IPFDonor_35kSortedNuclei_cDNA","UCinc_UMich_IPF0165_Lung_IPFDonor_14kDirectNuclei_cDNA","UCinc_UMich_IPF0165_Lung_IPFDonor_35kSortedNuclei_cDNA","UCinc_UNC_KKCFFT044N_Lung_HealthyDonor_14kDirectNuclei_cDNA","UCinc_UNC_KKCFFT044N_Lung_HealthyDonor_35kSortedNuclei_cDNA","UCinc_UFM_L0171_Lung_IPFDonor_FACSCounted15kNuclei_cDNA","UCinc_UFM_L0172_Lung_IPFDonor_FACSCounted15kNuclei_cDNA","UCinc_UNC_DD0025N_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA","UCinc_UNC_DD026P_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA","UCinc_UC_001_030_Lung_IPFDonor_NoSort15kNuclei_cDNA","UCinc_UFM_IPF0175_Lung_IPFDonor_NoSort15kNuclei_cDNA","UCinc_UFM_L0173_Lung_IPFDonor_NoSort15kNuclei_cDNA","UCinc_UFM_L201_Lung_IPFDonor_NoSort15kNuclei_cDNA","UCinc_UMich_L182_Lung_IPFDonor_NoSort15kNuclei_cDNA","UCinc_UNC_DD016N_Lung_HealthyDonor_NoSort15kNuclei_cDNA","UCinc_UNC_DD039N_Lung_HealthyDonor_NoSort15kNuclei_cDNA","UCinc_UNC_DD057N_Lung_HealthyDonor_NoSort15kNuclei_cDNA","UCinc_UC_001_007_Lung_IPFDonor_FACSCounted15kNuclei_cDNA","UCinc_UC_001_021_Lung_IPFDonor_FACSCounted15kNuclei_cDNA","UCinc_UC_001_026_Lung_IPFDonor_FACSCounted15kNuclei_cDNA","UCinc_UC_UC_1_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA","UCinc_UFM_L195_Lung_IPFDonor_FACSCounted15kNuclei_cDNA","UCinc_UNC_DD003O_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA","UCinc_UNC_DD020P_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA","UCinc_CCHMC_HL_CCHMC_01_Lung_HealthyDonor_FACSCounted16kNuclei_cDNA","UCinc_UFM_L210_Lung_IPFDonor_FACSCounted16kNuclei_cDNA","UCinc_UMich_UMILD146_LLL0212_Lung_IPFDonor_FACSCounted16kNuclei_cDNA","UCinc_UMILD0145_L0211_Lung_IPFDonor_FACSCounted16kNuclei_cDNA")

# To read all the seurat objects
for(i in 1:29){
  a <- readRDS(file = paste(flag[i], ".RDS", sep = ""))
  assign(paste(flag[i]), a)
}

# To save all the seurat objects
for(i in 1:29){
  a <- get(paste(flag[i]))
  saveRDS(a, file = paste(flag[i], ".RDS", sep = ""))
}

# Adding condition and renaming the cell ID's
for(i in 1:29){
  a <- get(paste(flag[i]))
  a <- RenameCells(a, add.cell.id = flag[i])
  a@meta.data$sample <- paste(flag[i])
  if(grepl("IPF", flag[i])){
    a@meta.data$condition <- "IPF"
  }
  else {
    a@meta.data$condition <- "healthy"
  }
  assign(paste(flag[i]), a)
}

# Calculating mito % and doing the general workflow for each object
for(i in 1:29){
  a <- get(paste(flag[i]))
  a[["percent.mt"]] <- PercentageFeatureSet(a, pattern = "^MT-")
  a <- subset(a, subset = nFeature_RNA > 1000 & percent.mt < 20)
  a <- NormalizeData(a, verbose = FALSE)
  a <- FindVariableFeatures(a, verbose = FALSE)
  a <- ScaleData(a)
  a <- RunPCA(a, dims = 1:30)
  assign(paste(flag[i]), a)
}

# Merging all objects
merged_seurat <- merge(UCinc_UFM_L169_Lung_IPFDonor_14kDirectNuclei_cDNA, y = c(UCinc_UFM_L169_Lung_IPFDonor_35kSortedNuclei_cDNA, UCinc_UMich_IPF0165_Lung_IPFDonor_14kDirectNuclei_cDNA, UCinc_UMich_IPF0165_Lung_IPFDonor_35kSortedNuclei_cDNA, UCinc_UFM_L0171_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UFM_L0172_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UC_001_030_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UFM_IPF0175_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UFM_L0173_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UC_001_007_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UFM_L201_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UMich_L182_Lung_IPFDonor_NoSort15kNuclei_cDNA, UCinc_UC_001_021_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UC_001_026_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UFM_L195_Lung_IPFDonor_FACSCounted15kNuclei_cDNA, UCinc_UFM_L210_Lung_IPFDonor_FACSCounted16kNuclei_cDNA, UCinc_UMich_UMILD146_LLL0212_Lung_IPFDonor_FACSCounted16kNuclei_cDNA, UCinc_UMILD0145_L0211_Lung_IPFDonor_FACSCounted16kNuclei_cDNA,UCinc_UNC_KKCFFT044N_Lung_HealthyDonor_14kDirectNuclei_cDNA, UCinc_UNC_KKCFFT044N_Lung_HealthyDonor_35kSortedNuclei_cDNA, UCinc_UNC_DD0025N_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_UNC_DD026P_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_UNC_DD016N_Lung_HealthyDonor_NoSort15kNuclei_cDNA, UCinc_UNC_DD039N_Lung_HealthyDonor_NoSort15kNuclei_cDNA, UCinc_UNC_DD057N_Lung_HealthyDonor_NoSort15kNuclei_cDNA, UCinc_UC_UC_1_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_UNC_DD003O_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_UNC_DD020P_Lung_HealthyDonor_FACSCounted15kNuclei_cDNA, UCinc_CCHMC_HL_CCHMC_01_Lung_HealthyDonor_FACSCounted16kNuclei_cDNA))

# General workflow for merged
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 1000 & percent.mt < 20)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 3000)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, npcs = 30)

# Harmony
merged_seurat_harmony <- merged_seurat %>% RunHarmony(group.by.vars = "sample", plot_convergence = F)
merged_seurat_harmony <- RunUMAP(merged_seurat_harmony, reduction = "harmony", dims = 1:20)
merged_seurat_harmony <- FindNeighbors(merged_seurat_harmony, reduction = "harmony", dims = 1:20)
merged_seurat_harmony <- FindClusters(merged_seurat_harmony, resolution = 0.6)

# Plots
FeaturePlot(combined.integrated, order = T, features = c("WT1"), cols = c("gray80", "red"), raster = FALSE, label = T, split.by = "condition")
DotPlot(combined.integrated, features = c("PLVAP","SELE","COL15A1","MPZL2","POSTN"))
VlnPlot(combined.integrated, features = c("SLC40A1"), split.by = "condition", raster = F, combine = F, pt.size = 0.5, sort = F, add.noise = T)
DimPlot(combined.integrated, reduction = "umap", label = T, raster = F)

# Heatmap has an issue, doesn't show any output
DoHeatmap(merged_seurat_harmony, features = c("WT1", "TP63", "SLC40A1", "KRT7"))


