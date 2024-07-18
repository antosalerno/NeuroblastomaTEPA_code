# Evaluate tumor heterogeneity #
## author: Antonietta Salerno
## date: 07/02/2023

library("Seurat") # Run with v4
library("MAST")
library("writexl")
library(openxlsx)
library(RColorBrewer)
library("EnhancedVolcano")
library("fgsea")
library(dplyr)
library(devtools)
library(msigdbr)
library(readr)
library(stringr)
library("org.Mm.eg.db", character.only = TRUE)
library(clusterProfiler)
library(ggplot2)
library(ggcharts)
library(glmGamPoi)

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
#tumor <- LoadSeuratRds("TEPA_results/S00_tumor.Rds") #13560 cells and 21388 genes

#### 1 - QC and filtering of tumor cells ####

png("TEPA_plots/S05_umapTumorNotInt.png", w = 2000, h = 2000, res = 300)
#pdf(qq("TEPA_final_figures/S05_umapTumorNotInt.pdf"), w = 6, h = 5)
DimPlot(object = tumor, pt.size = 0.000000005, group.by = "condition",
        reduction = 'umap', ncol = 1,
        label = FALSE, cols = cond_col) +
  ggtitle("Tumor cells in TEPA vs Control") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

tumor <- NormalizeData(tumor) 
# Add percent mito data in seurat seuset_tumorect
tumor[["percent.mt"]] <- PercentageFeatureSet(tumor, pattern = "^mt.")
# Add percent ribo genes
tumor[["percent.ribo"]] <- PercentageFeatureSet(tumor, pattern = "^Rp.")

# Visualize the distribution of mitochondrial gene expression detected per cell
png("TEPA_plots/S05_tumorPreFilterQC_mito.png", h = 3000, w = 4200, res = 300)
tumor@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)
dev.off()

png("TEPA_plots/S05_tumorPreFilterQC_ribo.png", h = 3000, w = 4200, res = 300)
tumor@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 15)
dev.off()

tumor <- subset(tumor, subset = nFeature_RNA > 200 &
                  percent.mt < 25 & percent.ribo < 15) # 13535 cells

png("TEPA_plots/S05_tumorPostFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(tumor, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
# a little difference in the distribution -> batch effect 
dev.off()

png("TEPA_plots/S05_tumorPostFilterQC_scatter.png", h = 3000, w = 4200, res = 300)
plot1 <- FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 # We can spot two separate data clouds because of the bimodal distribution of the data: control vs treated cells
dev.off()

### Batch effect correction ###

samples.list <- SplitObject(tumor, split.by = "condition")

# Normalize and identify variable features for each dataset independently (Treatment vs Control)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

for (i in 1:length(samples.list)) {
  samples.list[[i]] <- SCTransform(samples.list[[i]], vst.flavor = "v2", verbose = TRUE) 
}

features <- SelectIntegrationFeatures(object.list = samples.list)
tumor.anchors <- FindIntegrationAnchors(object.list = samples.list, anchor.features = features)
tumor.combined <- IntegrateData(anchorset = tumor.anchors, layer=)

#### 2 - Dimensionality reduction ####

# Scaling
seuset_tumor <- ScaleData(tumor.combined, verbose = FALSE)
seuset_tumor <- RunPCA(seuset_tumor, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_tumor@reductions$pca@stdev / sum(seuset_tumor@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 73.24 % of variability

seuset_tumor <- FindNeighbors( seuset_tumor, graph.name = "clust", dims = 1:60, reduction = 'pca')
seuset_tumor <- FindClusters(seuset_tumor, graph.name = "clust", resolution = 0.3, algorithm = 1) # original plot is 0.7
seuset_tumor <- RunUMAP(seuset_tumor, dims = 1:60, reduction = "pca", verbose = FALSE)

#png("TEPA_plots/S05_tumorUmapClust.png", w = 4000, h = 4000, res = 350)
pdf("TEPA_plots/S05_tumorUmapClust.pdf", w = 8, h = 8)
DimPlot(seuset_tumor, pt.size = 0.5, reduction = 'umap', ncol = 1, 
             group.by = c("seurat_clusters", "condition"), label = FALSE) +
  ggtitle(paste(as.character(nrow(seuset_tumor@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

#png("TEPA_plots/S05_tumorUmapClusters.png", w = 2000, h = 2000, res = 200)
pdf("TEPA_final_figures/S05_tumorUmapClusters_v5.pdf", w = 9, h = 7)
DimPlot(seuset_tumor, pt.size = 0.5, reduction = 'umap', ncol = 1, cols = tum_col,
        group.by = c("seurat_clusters"), label = TRUE) +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

png("TEPA_plots/S05_tumorUmapClusters.png", w = 2000, h = 2000, res = 200)
#pdf(qq("TEPA_final_figures/S05_tumorUmapCond.pdf"), w = 10, h = 7)
DimPlot(seuset_tumor, pt.size = 0.5, reduction = 'umap', ncol = 1, cols = cond_col,
        group.by = c("condition"), label = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

#SaveSeuratRds(seuset_tumor, "TEPA_results/S05_seusetTumorClu.Rds")

#### 3 - Clustering annotation ####

#seuset_tumor <- LoadSeuratRds("TEPA_results/S05_seusetTumorClu.Rds")


### 3.1 Inter-cluster DEA: get marker genes ###

seuset_tumor <- JoinLayers(seuset_tumor)

DefaultAssay(seuset_tumor) <- "integrated"

save = "S05_tumorMarkers_v5"
Idents(seuset_tumor) <- "seurat_clusters"


## A - Find markers for every cluster compared to all remaining cells
tumor.markers <- FindAllMarkers(seuset_tumor, 
                                only.pos = FALSE, 
                                min.pct = 0.5, 
                                min.diff.pct = 0.2,
                                logfc.threshold = 0.3, 
                                test.use="MAST",
                                latent.vars="orig.ident")

write.csv(tumor.markers, file=paste0("TEPA_results/", save,".csv"))

# Save results in different excel sheets 
clusters = levels(Idents(seuset_tumor))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = tumor.markers[tumor.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file=paste0("TEPA_results/",save,".xlsx"), overwrite = TRUE)

# B - Add module score to the Seurat object for each cluster ###
tumor.markers <- read.csv("TEPA_results/S05_tumorMarkers_v5.csv")

seuset_tumor <- createSets(markers = tumor.markers, 
                           obj = seuset_tumor, id = "seurat_clusters")

#png("TEPA_plots/S05_tumorClustersAnn.png", h = 3000, w = 4500, res = 300)
pdf("TEPA_final_figures/S05_tumorClustersAnn_v5.pdf", h = 7, w = 9)
clusters = unique(Idents(seuset_tumor))
patchwork::wrap_plots(FeaturePlot(seuset_tumor, ncol = 2, combine = TRUE, pt.size = 1, 
                                  features = as.character(clusters), label = TRUE, repel = TRUE)) & theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off() 

## C - Identify most important markers per cluster ###

markers <- getTopMarkers(tumor.markers, 5)

#png("TEPA_plots/S05_tumorDotPlot.png", h = 2000, w = 2500, res = 300)
pdf("TEPA_final_figures/S05_tumorDotPlot_v5.pdf", h = 5, w = 6)
DotPlot(seuset_tumor, features = unique(markers), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

## D - Plot Volcano of each cluster vs all the others: 
save = "S05_tumorMarkers"
plotVolcano(clusters, res = tumor.markers, type = "markers", immune = F, log2FC = 0.5, save = save)

#### 3.2 - Intra-cluster DEA with annotated dataset - Treatment vs Control ####
Idents(seuset_tumor) <- "seurat_clusters"
seuset_tumor$clusters.tepa <- paste(Idents(seuset_tumor), seuset_tumor$condition, sep = "_")
Idents(seuset_tumor) <- "clusters.tepa"
DefaultAssay(seuset_tumor) <- "RNA"

save = "S05_tumorCond_"
sheets <- list()
for (cluster in unique(seuset_tumor$seurat_clusters)){
  try({
    ident1 <- paste0(cluster,"_Treatment")
    ident2 <- paste0(cluster,"_Control")
    message(paste0("DEA control vs treatment for cluster: ", cluster))
    condition.diffgenes <- FindMarkers(seuset_tumor, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       #logfc.threshold = 0.25, 
                                       only.pos = FALSE, verbose = FALSE,
                                       #latent.vars="orig.ident",
                                       #min.cells.feature = 1, min.cells.group = 1, 
                                       test.use="MAST")
    condition.diffgenes$p_val_adj = p.adjust(condition.diffgenes$p_val, method='BH')
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    write.csv(condition.diffgenes, file=paste0("TEPA_results/", save,"DEA_",sub(" ", "_", cluster),".csv"))
  })
}
# Needed for manual curation
openxlsx::write.xlsx(sheets, paste0("TEPA_results/", save,"DEA.xlsx"), rowNames=TRUE)

# Plot Volcano DEA by condition

Idents(seuset_tumor) <- "seurat_clusters"
clusters = unique(Idents(seuset_tumor))

save = "S05_tumor_"
plotVolcano(clusters, log2FC = 0.5, save = save, immune = F)

#### 3.3 DEA Treatment vs Control bulk dataset ####
DefaultAssay(seuset_tumor) <- "RNA"
seuset_tumor <- JoinLayers(seuset_tumor)

counts <- GetAssayData(seuset_tumor, assay = "RNA", layer = "counts")
counts <- counts[-(which(rownames(counts) == "Xist")),]
seuset_tumor <- subset(seuset_tumor, features = rownames(counts))

Idents(seuset_tumor) <- "condition"

seuset_tumor <- PrepSCTFindMarkers(seuset_tumor)

res <- FindMarkers(seuset_tumor, 
                                ident.1 = "Treatment", ident.2 = "Control",
                                only.pos = FALSE, verbose = FALSE, layer="counts",
                                assay = "RNA",
                                latent.vars="orig.ident",
                                min.pct = 0.5, 
                                min.diff.pct = 0.2,
                                test.use="MAST")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/S05_tumorBulkDEA_MAST_v5.csv"))

log2FC = 0.5
save = "S05_tumorBulk_v5"
res <- res %>% 
  filter(p_val_adj < 0.05);
res[order(-res$avg_log2FC),]

markers_DEA <- rownames(rbind(head(res,20), tail(res,20)))
  
#png("TEPA_plots/S05_tumorBulk_DEA_FeaturePlot.png", h = 2000, w = 2500, res = 300)
pdf("TEPA_final_figures/S05_tumorBulk_DEA_FeaturePlot.pdf", h = 2, w = 10)
DotPlot(seuset_tumor, features = unique(markers_DEA), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()
# Stmn3 = Exhibits microtubule-destabilizing activity, which is antagonized by STAT3
# Zfp580 = cellular response to hydrogen peroxide; positive regulation of cell migration

p <- EnhancedVolcano(res, subtitle = "",
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3.5, 2),
                     title = "DEA Tumor TEPA vs Control",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 3,
                     col = c('lightgrey', 'pink', "#ADD8E6", 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     cutoffLineType = 'blank',
                     cutoffLineCol = 'black',
                     cutoffLineWidth = 0.8,
                     hline = c(10e-10, 10e-50, 10e-100, 10e-300),
                     hlineCol = c('grey0', 'grey25','grey50','grey75'),
                     hlineType = 'longdash',
                     hlineWidth = 0.8,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE,
                     widthConnectors = 0.3,
                     colConnectors = 'gray51', 
                     maxoverlapsConnectors = 170,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes')) + 
  coord_flip(expand = TRUE, clip = "on") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(0,400))
ggsave(p, file=paste0("TEPA_final_figures/", save, "DEA.pdf"), width = 28, height = 17, units = "cm")

#### 4 - Gene Set Enrichment Analysis ####

# Read all sets and store them in a list
sets <- lapply(paste0("TEPA_data/", c(
  "mh.all.v2022.1.Mm.symbols.gmt",
  "REACTOME_NEUTROPHIL_DEGRANULATION.v2022.1.Mm.gmt",
  "GOBP_CELL_MOTILITY.v2023.2.Mm.gmt"
  # "GOBP_NEUTROPHIL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_CHEMOTAXIS.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_DEGRANULATION.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_DIFFERENTIATION.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_EXTRAVASATION.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_MEDIATED_CYTOTOXICITY.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_MEDIATED_IMMUNITY.v2023.2.Mm.gmt",
  # "GOBP_NEUTROPHIL_MIGRATION.v2023.2.Mm.gmt",
  # "GOBP_CELLULAR_RESPONSE_TO_COPPER_ION.v2023.2.Mm.gmt",
  # "GOBP_COPPER_ION_HOMEOSTASIS.v2023.2.Mm.gmt",
  # "GOBP_COPPER_ION_IMPORT.v2023.2.Mm.gmt",
  # "GOBP_COPPER_ION_TRANSMEMBRANE_TRANSPORT.v2023.2.Mm.gmt",
  # "GOBP_COPPER_ION_TRANSPORT.v2023.2.Mm.gmt",
  # "GOBP_DETOXIFICATION_OF_COPPER_ION.v2023.2.Mm.gmt",
  # "GOBP_RESPONSE_TO_COPPER_ION.v2023.2.Mm.gmt",
  # "GOMF_COPPER_CHAPERONE_ACTIVITY.v2023.2.Mm.gmt",
  # "GOMF_COPPER_ION_BINDING.v2023.2.Mm.gmt",
  # "GOMF_COPPER_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY.v2023.2.Mm.gmt"
)), read.gmt)

# Create named list with ontology term and containing genes
sgsea <- fgsea_sets(sets)


Idents(seuset_tumor) <- "condition"
clusters = unique(Idents(seuset_tumor))
save = "S05_tumorEnrichment_Cond"
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, input = "tumor")

#### 5 - Look at expression of specific genes ####

# Search all isoforms of gene of interest
grep(pattern = "Ifi2", 
     x = rownames(x = seuset_tumor@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

#png("TEPA_plots/S05_tumorMt1.png", h = 2000, w = 2000, res = 200)
pdf("TEPA_final_figures/S05_tumorMt1.pdf", h = 4, w = 4)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  c("Mt1"),
        cols = cond_col, ncol = 1, pt.size = 0) + 
  ylim(-1,5) +
  geom_signif(xmin = 1, xmax = 2, y_position = 4, annotations="***")+
  geom_boxplot(width=.1, outlier.size = 0.5)
dev.off()

#png("TEPA_plots/S05_tumorMt2.png", h = 2000, w = 2000, res = 200)
pdf("TEPA_final_figures/S05_tumorMt2.pdf", h = 4, w = 4)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  c("Mt2"),
        cols = cond_col, ncol = 1, pt.size = 0) + 
  ylim(-1,5) +
  geom_signif(xmin = 1, xmax = 2, y_position = 4, annotations="***")+
  geom_boxplot(width=.1, outlier.size = 0.5)
dev.off()

#png("TEPA_plots/S05_tumorMycn.png", h = 2000, w = 2000, res = 200)
pdf(qq("TEPA_final_figures/S05_tumorMycn.pdf"), h = 4, w = 4)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  "Mycn", cols = cond_col,
        ncol = 1, pt.size = 0) + 
  ylim(0,6) +
  #geom_signif(xmin = 1, xmax = 2, y_position = 5.25, annotations="***") +
  geom_boxplot(width=.1, outlier.size = 0.5)
dev.off()

png("TEPA_plots/S05_tumorH1fx.png", h = 2000, w = 2000, res = 200)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  "H1fx", cols = cond_col,
        ncol = 1, pt.size = 0.000005) + 
  ylim(0,6) +
  geom_signif(xmin = 1, xmax = 2, y_position = 5.25, annotations="***")
dev.off()


