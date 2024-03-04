###### Subset tumour cells from immune cells ##
## author: Antonietta Salerno
## date: 21/12/2022

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

library("Seurat")
library("SeuratObject")
library("ggplot2")
library(dplyr)
library(RColorBrewer)
library(patchwork)

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")

#### 1 - Create Seurat object with all the samples ####

#d7c : day 7 - control
counts <- read.table("TEPA_data/Combined_D7_CONTROL_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
CF <- CreateSeuratObject(counts = t(counts), project="CF", min.cells = 3, min.features = 500)
#d7cm : day 7 - control myeloid 
counts <- read.table("TEPA_data/Combined_D7_CONTROLMYE_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
CM <- CreateSeuratObject(counts = t(counts), project="CM", min.cells = 3, min.features = 500) 
#d7t : day 7 - tepa
counts <- read.table("TEPA_data/Combined_D7_TEPA_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
TF <- CreateSeuratObject(counts = t(counts), project="TF", min.cells = 3, min.features = 500)
#d7Tm : day 7 - tepa myeloid 
counts <- read.table("TEPA_data/Combined_D7_TEPAMYE_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
TM <- CreateSeuratObject(counts = t(counts), project = "TM",min.cells = 3, min.features = 500)

# Create list object to be merged in a large Seurat object
seuset <- list(CM = CM, TF = TF, TM = TM)

combined <- merge(
  x = CF,
  y = seuset,
  add.cell.ids = c("CF",'CM','TF', 'TM')
)

seuset <- combined
seuset<-JoinLayers(seuset)
seuset$condition <- ifelse(test = seuset$orig.ident %in% c("CF", "CM"), yes = "Control", no = "Treatment")
seuset$sampleType <- ifelse(test = seuset$orig.ident %in% c("TM", "CM"), yes = "Myeloid", no = "Full")

# table(seuset$orig.ident)
SaveSeuratRds(seuset, "TEPA_results/S00_rawcounts.Rds")

#### 2 - Data pre-processing ####
seuset <- LoadSeuratRds("TEPA_results/S00_rawcounts.Rds")

png("TEPA_plots/S00_preFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(seuset, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S00_preFilterQC_Mycn_Cd45.png", h = 3000, w = 4200, res = 300)
FeatureScatter(seuset, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 0.0005, slot="counts")
dev.off()

# Log-normalisation
seuset <- NormalizeData(seuset) 

# We only carry out filtering for immune cells after separating the tumour

#### 3 - Dimensionality reduction ####

seuset <- FindVariableFeatures(seuset, 
                               selection.method = "vst",
                               nfeatures = 2000)
seuset <- ScaleData(seuset, verbose = FALSE)
seuset <- RunPCA(seuset, npcs = 100, verbose = FALSE) 

# Determine percent of variation associated with each PC
pct <- seuset@reductions$pca@stdev / sum(seuset@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=40) # Select 40 PCs to retain 62.76% of variability

seuset <- FindNeighbors(object = seuset, dims = 1:40, reduction = 'pca')
seuset <- RunUMAP(seuset, dims = 1:40, reduction = "pca", verbose = FALSE)

png("TEPA_plots/S00_umapExplore.png", w = 4000, h = 2000, res = 300)
#pdf("TEPA_final_figures/S00_umapExplore.pdf", w = 12, h = 6)
cbpal <- colorRampPalette(RColorBrewer::brewer.pal(4,"Paired"))(4)
DimPlot(object = seuset, pt.size = 0.05, reduction = 'umap', ncol = 2,
        group.by = c("orig.ident", "condition"), label = FALSE, cols = cbpal) +
  ggtitle(paste(as.character(nrow(seuset@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

png("TEPA_plots/S00_umapCounts.png", w = 4000, h = 2000, res = 300)
#pdf("TEPA_final_figures/S00_umapCounts.pdf", w = 12, h = 6)
FeaturePlot(seuset, features = c("nCount_RNA", "nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

png("TEPA_plots/S00_umapMycnPtprc.png", h = 2000, w = 4000, res = 300)
#pdf("TEPA_final_figures/S00_umapMycnPtprc.pdf", w = 12, h = 6)
FeaturePlot(seuset, ncol = 2, pt.size = 0.0005,
            features = c("Mycn", "Ptprc"), label = TRUE, repel = TRUE) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

SaveSeuratRds(seuset, "TEPA_results/S00_seusetRed.Rds")

#### 4 - Identify tumor cells ####

#seuset <- LoadSeuratRds("TEPA_results/S00_seusetRed.Rds")

table(seuset$orig.ident) # pre-subsetting: CF: 9100  CM: 7377  TF: 9412 TM: 11174 

png("TEPA_plots/S00_MycnDistribution.png", h = 3000, w = 4200, res = 300)
X <- colSums(seuset["Mycn",])
hist_dens(X, breaks = 100)
dev.off()

tumorCells1 <- subset(seuset[,seuset$sampleType == "Full"], Mycn > 0 & Ptprc == 0)
tumorCells2 <- subset(tumorCells1[,tumorCells1$sampleType == "Full"], 
                     Cd3d > 0 | Cd3e > 0 | Cd3g > 0 | Cd8a > 0 | Cd8b1 > 0, invert = TRUE) # 13560 cells: 169 lymphocytic cells rescued
immuneCells <- setdiff(Cells(tumorCells1), Cells(tumorCells2)) # rescued lymphocytes
seuset$class <- ifelse(test = colnames(seuset) %in% append(Cells(tumorCells1), Cells(tumorCells2)), "tumorCells", "noClass")

png("TEPA_plots/S00_tumorCells.png", h = 2000, w = 4000, res = 300)
FeaturePlot(seuset[,seuset$class == "tumorCells"], ncol = 2, pt.size = 0.005,
            features = c("Mycn","Ptprc"), label = TRUE, repel = TRUE) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

immuneCells <- append(immuneCells, Cells(subset(seuset[,seuset$sampleType == "Full"], Ptprc > 0 & Mycn == 0)))
immuneCells <- append(immuneCells, Cells(subset(seuset[,seuset$sampleType == "Myeloid"])))
seuset$class <- ifelse(test = colnames(seuset) %in% immuneCells, "immuneCells", seuset$class)
seuset$class <- ifelse(test = colnames(seuset) %in% Cells(append(tumorCells2,immuneCells)), "noClass", seuset$class)

pt <- table(Idents(seuset), seuset$class)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

png("TEPA_plots/S00_classMemberMycnPtprc.png", w=2500,h=2500, res=300)
#pdf("TEPA_final_figures/S00_classMemberIdent.pdf", w = 6, h = 6)
cbpal <- colorRampPalette(RColorBrewer::brewer.pal(4,"Paired"))(4)
pt1 <- DimPlot(object = seuset, pt.size = 0.05, reduction = 'umap', ncol = 1,
        group.by = c("orig.ident"), label = FALSE, cols = cbpal) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) 

pt2 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity")+
  xlab("Class") +
  ylab("Cell counts") +
  #ggtitle("Class membership of cell counts") +
  scale_fill_manual(values = brewer.pal(4, "Paired")) +
  theme(legend.title = element_blank())

(pt1 + pt2) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

dev.off()

png("TEPA_plots/S00_noClass_MycnPtprc.png", h = 3000, w = 4200, res = 300)
FeatureScatter(seuset[,seuset$class == "noClass"], feature1 = "Mycn", feature2 = "Ptprc", pt.size = 1)
dev.off()

SaveSeuratRds(seuset, "TEPA_results/S00_seusetClass.Rds")

### 5 -  Select only immune cells ####

#seuset <- LoadSeuratRds("TEPA_results/S00_seusetClass.Rds")

# remove noClass 
seuset <- subset(seuset, class != "noClass")

png("TEPA_plots/S00_PrnpByClass.png", h = 3000, w = 4200, res = 500)
VlnPlot(seuset, features = c("Prnp"), split.by = "condition" , group.by = "class",
        ncol = 1, pt.size = 0.0005)
dev.off()

png("TEPA_plots/S00_Mt1ByClass.png", h = 3000, w = 4200, res = 500)
VlnPlot(seuset, features = c("Mt1"), split.by = "condition" , group.by = "class",
        ncol = 1, pt.size = 0.0005) +
  ylim(0,5) 
dev.off()

# Tumor cells
tumor <- seuset[,seuset$class == "tumorCells"]
table(tumor$orig.ident) 
SaveSeuratRds(tumor, "TEPA_results/S00_tumor.Rds")

# Immune cells
immune <- seuset[,seuset$class == "immuneCells"]
table(immune$orig.ident) 
# post-subsetting all samples: CF: 2147 CM: 7377  TF: 1484 TM: 11174

png("TEPA_plots/S00_immuneCells_MycnPtprc.png", h = 3000, w = 4200, res = 300)
FeatureScatter(immune, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 1)
dev.off()

png("TEPA_plots/S00_postFilterImmune_violin.png", h = 3000, w = 4200, res = 300)
#pdf("TEPA_final_figures/S00_postFilterImmune_violin.pdf", w = 15, h = 8)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), 
        ncol = 3, pt.size = 0.00000000005) &
  scale_fill_manual(values = brewer.pal(4, "Paired"))
dev.off()

SaveSeuratRds(immune,"TEPA_results/S00_immune.Rds")


## unload currently loaded packages for next session
#lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))


