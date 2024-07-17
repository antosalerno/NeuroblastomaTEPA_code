###### Explore the cells not classified as tumor nor noclass ##
## author: Antonietta Salerno
## date: 04/03/2024

# In the full samples CF and TF there are some cells that are not yet classified, they express both Mycn and Ptprc

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
library("writexl")
library(openxlsx)
library('limma')
library("MAST")


setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")

seuset <- LoadSeuratRds("TEPA_results/S00_seusetClass.Rds")

noclass <- subset(seuset, class == "noClass")

# CF = 840, TF = 481

png("TEPA_plots/S00_postFilterNoClass_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(noclass, features = c("nFeature_RNA", "nCount_RNA", "Mycn", "Ptprc", "Vegfa", "Fap"), 
        ncol = 2, pt.size = 0.00000000005) &
  scale_fill_manual(values = brewer.pal(4, "Paired"))
dev.off()

# The treated samples have higher transcripts counts, more genes expressed but similar Mycn and Ptprc
# more Vegfa in treatment
# Mesenchymal stromal cells? Cancer Associated Fibroblasts?

noclass <- NormalizeData(noclass) 

# Add percent mito data in seurat object
noclass[["percent.mt"]] <- PercentageFeatureSet(noclass, pattern = "^mt.")
# Add percent ribo genes
noclass[["percent.ribo"]] <- PercentageFeatureSet(noclass, pattern = "^Rp.")

png("TEPA_plots/S07_preFilterQC_scatter.png", w = 6000, h = 2000, res = 300)
plot1 <- FeatureScatter(noclass, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(noclass, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.0005)
plot3 <- FeatureScatter(noclass, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 + plot3
dev.off()

# Visualize the distribution of mitochondrial and ribosomal gene expression detected per cell
png("TEPA_plots/S07_preFilterQC_mito.png", h = 3000, w = 4200, res = 300)
noclass@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)
dev.off()

png("TEPA_plots/S07_preFilterQC_ribo.png", h = 3000, w = 4200, res = 300)
noclass@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 30)
dev.off()

noclass <- subset(noclass, subset = percent.ribo < 10 & nCount_RNA > 200) # 1321 >> 672
#try first without dicarding cells with high percentage of mitochondrial genes (most of them)

png("TEPA_plots/S07_postFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(noclass, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S07_postFilterQC_Mycn_Cd45.png", h = 3000, w = 4200, res = 300)
FeatureScatter(noclass, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 0.0005)
dev.off()

noclass <- subset(noclass, subset = Mycn > 1.5 & Ptprc < 2.5) # 672 > 390 cells because of artifacts

png("TEPA_plots/S07_postFilterQC_Mycn_Cd45_filt.png", h = 3000, w = 4200, res = 600)
FeatureScatter(noclass, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 0.0005)
dev.off()

png("TEPA_plots/S07_postFilterQC_scatter.png", h = 3000, w = 4200, res = 300)
plot1 <- FeatureScatter(noclass, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(noclass, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.0005)
plot3 <- FeatureScatter(noclass, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 + plot3  # We can spot two separate data clouds because of the bimodal distribution of the data: tumour vs noclass cells
dev.off()

SaveSeuratRds(noclass,"TEPA_results/S07_noclassFilt.Rds")

#### 2 - Batch effect correction ####

#noclass <- LoadSeuratRds("S01_noclassFilt.Rds")

samples.list <- SplitObject(noclass, split.by = "condition")

# Normalize and identify variable features for each dataset independently (Treatment vs Control)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = samples.list)
noclass.anchors <- FindIntegrationAnchors(object.list = samples.list,anchor.features = features)
noclass.combined <- IntegrateData(anchorset = noclass.anchors)

SaveSeuratRds(noclass.combined, "TEPA_results/S07_noclassInt.Rds")


### 2 - Clustering ####

# Scaling
noclass <- ScaleData(noclass.combined, verbose = FALSE)
noclass <- RunPCA(noclass, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- noclass@reductions$pca@stdev / sum(noclass@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 68.64% of variability

noclass <- FindNeighbors(object = noclass, graph.name = "clust", dims = 1:60, reduction = 'pca')

noclass <- FindClusters(object = noclass, graph.name = "clust", resolution = 0.3) 
noclass <- RunUMAP(noclass, dims = 1:60, reduction = "pca", verbose = FALSE)

png("TEPA_plots/S07_umapExplore.png", w = 6000, h = 2000, res = 300)
#pdf(qq("TEPA_final_figures/S07_umapExplore.pdf"), w = 15, h = 8)
DimPlot(object = noclass, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("orig.ident", "condition", "seurat_clusters"), label = FALSE) +
  ggtitle(paste(as.character(nrow(noclass@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("TEPA_plots/S07_umapCounts.png", w = 4000, h = 2000, res = 300)
FeaturePlot(noclass, features = c("nCount_RNA", "nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

png("TEPA_plots/S07_umapClust.png", w = 4000, h = 2000, res = 300)
DimPlot(object = noclass, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("seurat_clusters"), split.by= "condition",label = TRUE) +
  ggtitle(paste(as.character(nrow(noclass@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

SaveSeuratRds(noclass, "TEPA_results/S07_noclassAnn.SeuratRds")


#### 3 - Inter-cluster DEA: get marker genes ####

seuset_noclass <- LoadSeuratRds("TEPA_results/S07_noclassAnn.SeuratRds")


DefaultAssay(seuset_noclass) <- "RNA"
Idents(seuset_noclass) <- "seurat_clusters"
clusters = levels(Idents(seuset_noclass))
seuset_noclass <- JoinLayers(seuset_noclass)

# A - Find markers for every cluster compared to all remaining cells
noclass.markers <- FindAllMarkers(seuset_noclass, 
                                 only.pos = FALSE, 
                                 min.pct = 0.5, 
                                 min.diff.pct = 0.2,
                                 logfc.threshold = 0.5, 
                                 test.use="MAST",
                                 latent.vars="orig.ident")
noclass.markers$p_val_adj = p.adjust(noclass.markers$p_val, method='BH')
write.csv(noclass.markers, "TEPA_results/S07_DEA_clusterMarkers.csv")

# Save results in different excel sheets 
clusters = levels(Idents(seuset_noclass))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = noclass.markers[noclass.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/S07_DEA_clusterMarkers.xlsx", overwrite = TRUE)

# B - Add module score to the Seurat object for each cluster ###

seuset_noclass <- createSets()

Idents(seuset_noclass) <- "seurat_clusters"
png("TEPA_plots/S07_noclassClustersAnnot.png", h = 3000, w = 4500, res = 200)
#pdf(qq("TEPA_final_figures/S07_noclassClustersAnnot.pdf"), h = 10, w = 14)
patchwork::wrap_plots(FeaturePlot(seuset_noclass, ncol = 5, combine = TRUE,
                                  features = as.character(clusters), label = FALSE, repel = TRUE)) &
  theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off() 

# C - Identify most important markers per cluster ###

markers <- getTopMarkers(noclass.markers, 5)

png("TEPA_plots/S07_noclassDotPlot.png", h = 2000, w = 2500, res = 300)
#pdf(qq("TEPA_final_figures/S07_noclassDotPlot.pdf"), h = 4, w = 11)
DotPlot(object = seuset_noclass, features = unique(markers), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 4, cols = c("blue","red")) + RotatedAxis() + #scale_x_reverse() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=8))
dev.off()

png("TEPA_plots/S07_noclassMarkersHeatmap.png", h = 5000, w = 6000, res = 300)
top10 <- noclass.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(object = subset(seuset_noclass, downsample = 500), size = 6, 
          assay = "integrated", features = top10$gene) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()



