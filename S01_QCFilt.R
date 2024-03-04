###### Dimensionality reduction and clustering of immune cells ##
## author: Antonietta Salerno
## date: 22/12/2022

library("Seurat")
library("ggplot2")

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
immune <- LoadSeuratRds("TEPA_results/S00_immune.Rds")

#### 1 - QC and filtering of immune cells ####

immune <- NormalizeData(immune) 

# Add percent mito data in seurat object
immune[["percent.mt"]] <- PercentageFeatureSet(immune, pattern = "^mt.")
# Add percent ribo genes
immune[["percent.ribo"]] <- PercentageFeatureSet(immune, pattern = "^Rp.")

png("TEPA_plots/S01_preFilterQC_scatter.png", w = 6000, h = 2000, res = 300)
plot1 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.0005)
plot3 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 + plot3
dev.off()

# Visualize the distribution of mitochondrial and ribosomal gene expression detected per cell
png("TEPA_plots/S01_preFilterQC_mito.png", h = 3000, w = 4200, res = 300)
immune@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)
dev.off()

png("TEPA_plots/S01_preFilterQC_ribo.png", h = 3000, w = 4200, res = 300)
immune@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 30)
dev.off()

immune <- subset(immune, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 &
                   percent.mt < 25 & percent.ribo < 20 & nCount_RNA < 3000) # 12127 cells left

png("TEPA_plots/S01_postFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S01_postFilterQC_Mycn_Cd45.png", h = 3000, w = 4200, res = 300)
FeatureScatter(immune, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 0.0005)
dev.off()

png("TEPA_plots/S01_postFilterQC_scatter.png", h = 3000, w = 4200, res = 300)
plot1 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.0005)
plot3 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 + plot3  # We can spot two separate data clouds because of the bimodal distribution of the data: tumour vs immune cells
dev.off()

SaveSeuratRds(immune,"TEPA_results/S01_immuneFilt.Rds", overwrite = TRUE)

#### 2 - Batch effect correction ####

immune <- LoadSeuratRds("S01_immuneFilt.Rds")

samples.list <- SplitObject(immune, split.by = "condition")

# Normalize and identify variable features for each dataset independently (Treatment vs Control)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = samples.list)
immune.anchors <- FindIntegrationAnchors(object.list = samples.list,anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

SaveSeuratRds(immune.combined, "TEPA_results/S01_immuneInt.Rds")





