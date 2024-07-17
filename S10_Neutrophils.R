##### Neutrophils sub-populations ####
## author: Antonietta Salerno
## date: 12/05/2024

library("Seurat")
library("ggplot2")
library("writexl")
library(openxlsx)
library(HGNChelper)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(leidenAlg)

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadSeuratRds("TEPA_results/S03_immuneDiff.Rds")


# Sub-clustering
neutro <- seuset_immune[,seuset_immune$celltypes == "Neutrophils"]
DefaultAssay(neutro) <- "integrated"

# Scaling
neutro <- ScaleData(neutro, verbose = FALSE, assay = 'integrated')
neutro <- RunPCA(neutro, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_immune@reductions$pca@stdev / sum(seuset_immune@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 70.26% of variability


neutro <- FindClusters(object = neutro, graph.name = "clust", resolution = 0.2, algorithm = 1)  #0.2 creates 2 sets and 0.4 creates 5


neutro <- RunUMAP(neutro, dims = 1:60, reduction = "pca", verbose = FALSE, assay = 'integrated')
DimPlot(neutro, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(neutro) <- "seurat_clusters"

#png("TEPA_plots/S10_umapExplore.png", w = 6000, h = 2000, res = 300)
pdf("TEPA_plots/S10_umapExplore.pdf", w = 8, h = 8)
DimPlot(object = neutro, reduction = 'umap', ncol = 1,
        group.by = c("seurat_clusters"), label = TRUE, label.size = 3, pt.size = 0.5) +
  ggtitle(paste(as.character(nrow(neutro@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# A - Find markers for every cluster compared to all remaining cells
Idents(neutro) <-"seurat_clusters"
neutro.markers <- FindAllMarkers(neutro, 
                                 only.pos = FALSE, 
                                 min.pct = 0.5, 
                                 min.diff.pct = 0.2,
                                 logfc.threshold = 0.25, 
                                 test.use="MAST",
                                 latent.vars=c("orig.ident", "condition"))
neutro.markers$p_val_adj = p.adjust(neutro.markers$p_val, method='BH')
write.csv(neutro.markers, "TEPA_results/S10_DEA_clusterMarkersNeutro.csv")

# Save results in different excel sheets 
clusters = levels(Idents(neutro))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = neutro.markers[neutro.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/S10_DEA_clusterMarkersNeutro.xlsx", overwrite = TRUE)

markers <- getTopMarkers(neutro.markers, 10)

#png("TEPA_plots/S10_NeutroDotPlot.png", h = 2000, w = 2500, res = 300)
pdf("TEPA_plots/S10_neutroDotPlot3.pdf", h = 4, w = 11)
DotPlot(object = neutro, features = unique(markers), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 4, cols = c("blue","red")) + RotatedAxis() + #scale_x_reverse() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=8))
dev.off()


#### Annotation ####

neutro$seurat_clusters <- as.factor(ifelse(test=neutro$seurat_clusters== "1", yes="N1", no='N2'))
Idents(neutro) <- neutro$seurat_clusters
SaveSeuratRds(neutro,"TEPA_results/S10_neutroAnn.Rds")

### Check if N1 more TEPA/Control 

neutro <- LoadSeuratRds("TEPA_results/S10_neutroAnn.Rds")


pdf("TEPA_plots/S10_N1_condition.pdf", w = 15, h = 8)
DimPlot(object = neutro, reduction = 'umap', ncol = 2,
        group.by = c("seurat_clusters", "condition"), label = TRUE, 
        label.size = 3, pt.size = 0.5) +
  ggtitle(paste(as.character(nrow(neutro@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("TEPA_plots/S10_Neutrophils_genes.pdf", h = 15, w = 15)
features=c("Mmp9", "S100a8", "Retnlg")
FeaturePlot(neutro, features = features, split.by = "condition")
dev.off()


### Compare conditions ####

neutro$neutro.tepa <- paste(Idents(neutro), neutro$condition, sep = "_")
#neutro$seurat_clusters <- Idents(neutro)
Idents(neutro) <- "neutro.tepa"
DefaultAssay(neutro) <- "RNA"

save = "S10_neutroCond_"
sheets <- list()
for (cluster in unique(neutro$seurat_clusters)){
  try({
    ident1 <- paste0(cluster,"_Treatment")
    ident2 <- paste0(cluster,"_Control")
    condition.diffgenes <- FindMarkers(neutro, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       logfc.threshold = 0.2, 
                                       min.pct = 0.5, 
                                       min.diff.pct = 0.2,
                                       only.pos = FALSE, verbose = FALSE,
                                       latent.vars="orig.ident",
                                       #min.cells.feature = 1, min.cells.group = 1, 
                                       test.use="MAST")
    condition.diffgenes$p_val_adj = p.adjust(condition.diffgenes$p_val, method='BH')
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    file=paste0("TEPA_results/", save, "DEA_", gsub(" |/", "_", cluster),".csv")
    write.csv(condition.diffgenes, file=file)
  })
}

# Needed for manual curation
openxlsx::write.xlsx(sheets, file = paste0("TEPA_results/", save, "DEA.xlsx"), rowNames=TRUE)

# Plot Volcano DEA by condition

Idents(neutro) <- "seurat_clusters"
clusters = unique(Idents(neutro))

plotVolcano(clusters, log2FC = 1, pval = 0.05, save = save, file = file)

# Plot N1/N2 signature

N1 <- c("Cebpb", "S100a8","S100a9","Ccl3", "Ifitm2", "Ifitm3", "Ifitm6", 
        "Acod1","Sell", "Prkcd","Retnlg","Mmp8", "Hif1a",  
        "Tnf", "Myd88", "Fas", "Cxcl3", "Isg15", "Isg20", "Arg2", 
        "Il15")

save = "S06_complexHeat_N1_FC_Neutro2"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 5)
sign_avgLogFCHeatMap(seuset_full, c(N1), immune = FALSE, celltype = "Neutrophils",
                     cluster = F, k = 0, legend = TRUE) 
dev.off()

N2 <- c("Tgfb1", "Ccl4", "Ccl5", "Il1r1", "Chil3", "Retnla", "Il1b","Il6ra", "Mmp9")

save = "S06_complexHeat_N2_FC_Neutro2"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 5)
sign_avgLogFCHeatMap(seuset_full, c(N2), immune = FALSE, celltype = "Neutrophils",
                     cluster = F, k = 0, legend = TRUE) 
dev.off()

pdf("TEPA_plots/S06_Ftl1_VlnPlot_Neutro.pdf", h = 6, w = 10)
Idents(neutro) <- "celltypes"
VlnPlot(neutro, features = "Ftl1", split.by = "condition", layer = "data", 
        ncol = 2, pt.size = 0.000005)
dev.off()

immune_migr <- c("Trem1", "Adam8", "Clec4e", "Cd14", "Ccr7", "Cxcr4", "Fpr1")

immune_extra <-c("Itgam", "Itgb2", "Cd177", "Itgal", "Cd44", "Itga2", "Itgb1", "Myo1f")

save = "S06_complexHeat_ExtraMigr_FC_Neutro"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 5)
sign_avgLogFCHeatMap(neutro, c(immune_extra, immune_migr), immune = TRUE, 
                     celltype = "Neutrophils",
                     cluster = TRUE, k = 1, legend = TRUE) 
dev.off()



