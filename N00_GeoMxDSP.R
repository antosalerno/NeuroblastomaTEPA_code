#### Analysis GeoMx DSP data ####
## author: Antonietta Salerno
## date: 07/03/2022
BiocManager::install("NanoStringNCTools")
BiocManager::install("GeomxTools")
BiocManager::install("GeoMxWorkflows")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(knitr)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggcharts)
library(MAST)
library(openxlsx)
library(patchwork)
library(clusterProfiler)
library("tibble")
library("Seurat")
library("SeuratDisk")
library("fgsea")
library("stringr")

### 1 - Load data ####
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")

datadir <- "TEPA_data"
DCCFiles <- dir(file.path(datadir, "DCC"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <- file.path(datadir, "annotationFull.xlsx")

data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template",
                               phenoDataDccColName = "Sample_ID",
                               protocolDataColNames = c("area", "roi", "Condition", 
                                                        "Core", "Infiltration",
                                                        "AOISurfaceArea", "AOINucleiCount",
                                                        "ROICoordinateX", "ROICoordinateY"),
                               experimentDataColNames = c("Core", "aoi"))

# Clean the column with Core (C,C3,C7,T3,T7)
pData(protocolData(data))$Core <- 
  str_split(pData(protocolData(data))$Core, "\\-", simplify=T)[,1]

# Add a column with Group (Control, T3, T7)
group <- pData(protocolData(data))$Core
control <-  group %in% c("C","C3","C7")
group[control] <- "C"
pData(protocolData(data))$Group <- group

### 2 - Study Design ####

count_mat <- dplyr::count(pData(protocolData(data)), Core, Group, Condition, Infiltration) %>%
  mutate(Condition = as.character(Condition)) %>% 
  mutate(Infiltration = as.character(Infiltration))

# gather the data and plot in order: 
test_gr <- gather_set_data(count_mat, 1:3)
test_gr$x <- factor(test_gr$x)
levels(test_gr$x) = c("Core","Group","Condition", "Infiltration")

# plot Sankey
save = "N00_sankeyCore"
p <- ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Infiltration), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "#E3B9B1", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") +
  annotate(geom = "segment", x = 3.25, xend = 3.25,
           y = 0, yend = 105, lwd = 2) +
  annotate(geom = "text", x = 3.19, y = 63, angle = 90, size = 5,
           hjust = 0.5, label = "100 ROIs")
#ggsave(p, file=paste0("TEPA_plots/", save, ".png"), width = 30, height = 30, units = "cm")
ggsave(p, file=paste0("TEPA_final_figures/", save, ".pdf"), width = 30, height = 30, units = "cm")

#### 3 - QC and preprocessing ####

### 3.1 - Shift counts to one###
data <- shiftCountsOne(data, useDALogic = TRUE)

### 3.2 - Flag low quality ROIs ###
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000) 
       percentTrimmed = 85,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%) ->75
       percentSaturation = 65, # Minimum sequencing saturation (50%)
       minNegativeCount = 1   # Minimum negative control counts (10) -> 1
  )   

data <- setSegmentQCFlags(data, qcCutoffs = QC_params)

### 3.3 - Flag low quality probes ###
data <- setBioProbeQCFlags(data, qcCutoffs = 
                             list(minProbeRatio = 0.1,percentFailGrubbs = 20), 
                           removeLocalOutliers = TRUE)

ProbeQCResults <- fData(data)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df

### 3.4 - Remove low quality ROIs and probes ###

passedQC <- 
  subset(data, 
         fData(data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(passedQC)

data <- passedQC 

### 3.5 - Create gene-level count data ###
# Objects must be aggregated to Target level data before coercing. 
# This changes the row (gene) information to be the gene name rather than the probe ID.

target_data <- aggregateCounts(passedQC)

# ### 3.6 - Limit of Quantification (Detection Rates) ###
# 
# # Define LOQ SD threshold and minimum value
# cutoff <- 2
# minLOQ <- 2
# 
# # Calculate LOQ per module tested 
# LOQ <- data.frame(row.names = colnames(target_data))
# LOQ <- pmax(minLOQ, pData(target_data)[, "NegGeoMean_Mm_R_NGS_WTA_v1.0"]
#             * pData(target_data)[, "NegGeoSD_Mm_R_NGS_WTA_v1.0"] ^ cutoff)
# 
# pData(target_data)$LOQ <- LOQ
# 
# ### - Filtering out either segments and/or genes with abnormally low signal ###
# LOQ_Mat <- c()
# ind <- fData(target_data)$Module == "Mm_R_NGS_WTA_v1.0"
# Mat_i <- t(esApply(target_data[ind, ], MARGIN = 1,
#                      FUN = function(x) {
#                        x > LOQ
#                      }))
# LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
# 
# # ensure ordering since this is stored outside of the geomxSet
# LOQ_Mat <- LOQ_Mat[fData(target_data)$TargetName, ]
# 
# ### By segment gene detection rates ###
# 
# # Save detection rate information to pheno data
# pData(target_data)$GenesDetected <- 
#   colSums(LOQ_Mat, na.rm = TRUE)
# pData(target_data)$GeneDetectionRate <-
#   pData(target_data)$GenesDetected / nrow(target_data)
# 
# target_data <-
#   target_data[, pData(target_data)$GeneDetectionRate >= .1]
# 
# dim(target_data)

### 3.6 - Normalisation 
### 3.6.1.  3rd quantile ###
norm_target_data <- normalize(target_data, norm_method="quant",
                                  desiredQuantile = .75, toElt = "q_norm")

# Try other normalisation methods: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9800292/

#### 4 - Coercion to Seurat object ####
seuset_nano <- as.Seurat(norm_target_data, normData = "q_norm", ident = "Group", 
                         coordinates = c("ROICoordinateX", "ROICoordinateY"))
Idents(seuset_nano) <- factor(x = Idents(seuset_nano), levels = c("C", "T3", "T7"))

# Remove T3 samples
seuset_nano <- subset(seuset_nano, Group != "T3")

seuset_nano@misc <- list()
SaveH5Seurat(seuset_nano, filename = "TEPA_results/N00_seusetNano.h5Seurat", overwrite = TRUE)

#head(seuset_nano@misc$QCMetrics$QCFlags) 

png("TEPA_plots/N00_countsROIs.png", h = 3000, w = 2500, res = 300)
VlnPlot(seuset_nano, features = "nCount_GeoMx", split.by = "Infiltration",
        pt.size = 0.1)
dev.off() # the number of genes is instead the same for all ROIs

png("TEPA_plots/N00_nucleiROIs.png", h = 3000, w = 2500, res = 300)
VlnPlot(seuset_nano, features = "AOINucleiCount", split.by = "Group",
        pt.size = 0.1)
dev.off() 

# ### Correct for batch effect ###
# 
# samples.list <- SplitObject(seuset_nano, split.by = "Group")
# 
# # Normalize and identify variable features for each dataset independently (Treatment vs Control)
# samples.list <- lapply(X = samples.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# features <- SelectIntegrationFeatures(object.list = samples.list)
# nano.anchors <- FindIntegrationAnchors(object.list = samples.list,
#                                        dims = 1:4, 
#                                        k.anchor = 1, k.filter = 1,
#                                        anchor.features = features)
# nano.combined <- IntegrateData(anchorset = nano.anchors)


# Then go to a new file and first do the tests (N01) and 

#### 5 - Dimensionality reduction ####
seuset_nano <- LoadH5Seurat("TEPA_results/N00_seusetNano.h5Seurat")

seuset_nano <- FindVariableFeatures(seuset_nano)
seuset_nano <- ScaleData(seuset_nano)
seuset_nano <- RunPCA(seuset_nano, assay = "GeoMx", verbose = FALSE, approx=FALSE)

# Determine percent of variation associated with each PC
pct <- seuset_nano@reductions$pca@stdev / sum(seuset_nano@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=50) # Select 50 PCs to retain 86.55% of variability

seuset_nano <- FindNeighbors(seuset_nano, reduction = "pca", dims = 1:50)
seuset_nano <- FindClusters(seuset_nano, resolution = 0.8)
seuset_nano <- RunUMAP(seuset_nano, reduction = "pca", dims = 1:50)

png("TEPA_plots/N00_clustROIs.png", h = 3000, w = 2500, res = 300)
DimPlot(seuset_nano, reduction = "umap", pt.size = 5,
        label = F, group.by = "seurat_clusters")
dev.off()
# 2 groups, would they reflect Core?

png("TEPA_plots/N00_umapExplore.png", w = 6000, h = 4000, res = 300)
DimPlot(object = seuset_nano, pt.size = 5, reduction = 'umap', ncol = 2, 
        group.by = c("Group", "Infiltration","Condition","seurat_clusters"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_nano@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off() # It looks like there's no pattern

#### 6 - Differential expression analysis ####
# Check in parallel expression of DEA genes in the different cell types
seuset_immune <- LoadH5Seurat("TEPA_results/S04_immuneDiff.h5Seurat")
seuset_nano@assays$GeoMx@scale.data <- scale(seuset_nano@assays$GeoMx@counts)

# 5.2 Infiltration vs Non-Infiltration given Treatment

Idents(seuset_nano) <- "Infiltration"
seuset_nanoTEPA <- subset(seuset_nano, Condition == "Treatment")

seuset_nanoTEPA@assays$GeoMx@scale.data <- scale(seuset_nanoTEPA@assays$GeoMx@counts)
res <- FindMarkers(seuset_nanoTEPA, ident.1 = "T", ident.2 = "F",
                   only.pos = FALSE, verbose = FALSE, assay= "GeoMx",
                   test.use="negbinom")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoInf_gCond_DEA.csv"))

log2FC = 0.5
save = "N00_nanoInf_gCond_"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3.5, 6),
                     ylim = c(0, 9),
                     title = "DEA Infiltrated vs Non-Infiltrated TME (TEPA-treated ROIs)",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 4,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,
                     colConnectors = 'gray51', maxoverlapsConnectors = 40,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
#ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")
ggsave(p, file=paste0("TEPA_final_figures/", save, "DEA.pdf"), width = 30, height = 25, units = "cm")

# ZNRF1 involved in ubiquitination of EGFR
# Mt2 down

ranked.genes<- res %>%
  filter(p_val_adj < 0.1) %>%
  arrange(desc(avg_log2FC))
ranked.genes <- rownames(ranked.genes)[!is.na(rownames(ranked.genes))]

png(paste0("TEPA_plots/",save,"DotPlot.png"), h = 2000, w = 2500, res = 300)
Idents(seuset_immune) <- "celltypes"
DotPlot(object = seuset_immune, features = ranked.genes[1:30], 
        scale=TRUE, dot.scale = 5,
        assay = "RNA", cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

png(paste0("TEPA_plots/",save,"Heatmap.png"), h = 4000, w = 6000, res = 300)
Idents(seuset_nano) <- "Infiltration"
DoHeatmap(object = subset(seuset_nanoTEPA, downsample = 500), size = 6, 
          group.by = "Infiltration",
          features = head(ranked.genes, 30)) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend()
dev.off()


# 5.4 Treatment vs Control given Infiltration

Idents(seuset_nano) <- "Condition"
seuset_nanoINF <- subset(seuset_nano, Infiltration == "T")
res <- FindMarkers(seuset_nanoINF, ident.1 = "Treatment", ident.2 = "Control",
                   only.pos = FALSE, verbose = FALSE,
                   #latent.vars="Core",
                   test.use="negbinom")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoCond_gInf_DEA.csv"))

log2FC = 0.5
save = "N00_nanoCond_gInf"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     ylim = c(0, 10),
                     title = "DEA Treatment vs Control (Infiltrated ROIs)",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 4,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,
                     colConnectors = 'gray51', maxoverlapsConnectors = 80,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
#ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")
ggsave(p, file=paste0("TEPA_final_figures/", save, "DEA.pdf"), width = 30, height = 25, units = "cm")

ranked.genes<- res %>%
  filter(p_val_adj < 0.1) %>%
  arrange(desc(avg_log2FC))
ranked.genes <- rownames(ranked.genes)[!is.na(rownames(ranked.genes))]

png(paste0("TEPA_plots/",save,"DotPlot.png"), h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_immune, features = ranked.genes[1:30], 
        scale=TRUE, dot.scale = 5, 
        assay = "RNA", cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

png(paste0("TEPA_plots/",save,"Heatmap.png"), h = 4000, w = 6000, res = 300)
Idents(seuset_nanoINF) <- "Condition"
DoHeatmap(object = subset(seuset_nanoINF, downsample = 500), size = 6, 
          features = head(ranked.genes, 40)) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

# 5.5 Treatment vs Control

Idents(seuset_nano) <- "Condition"
res <- FindMarkers(seuset_nano, ident.1 = "Treatment", ident.2 = "Control",
                   only.pos = FALSE, verbose = FALSE,
                   #latent.vars="Core",
                   test.use="negbinom")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoCond_DEA.csv"))

log2FC = 0.5
save = "N00_nanoCond"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     title = "DEA Treatment vs Control",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 4,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")

ranked.genes<- res %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > abs(0.3))
ranked.genes <- rownames(ranked.genes)[!is.na(rownames(ranked.genes))]

png(paste0("TEPA_plots/",save,"DotPlot.png"), h = 2000, w = 2500, res = 300)
Idents(seuset_immune) <- "celltypes"
DotPlot(object = seuset_immune, features = ranked.genes[1:30],  
        scale=TRUE, dot.scale = 5,
        assay = "RNA", cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

png(paste0("TEPA_plots/",save,"Heatmap.png"), h = 4000, w = 6000, res = 300)
DoHeatmap(object = subset(seuset_nano, downsample = 500), size = 6, 
          features = head(ranked.genes, 20)) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

SaveH5Seurat(seuset_nano, filename = "TEPA_results/N00_seusetNanoRed.h5Seurat", overwrite = TRUE)

#### 7 - Pathway enrichment Analysis ####
seuset_nano <- LoadH5Seurat("TEPA_results/N00_seusetNanoRed.h5Seurat")

### 7.1 - Cluster markers
### 7.2 - Infiltration vs Non-Infiltration given Treatment

sets1 <- read.gmt("TEPA_data/mh.all.v2022.1.Mm.symbols.gmt") # Mouse hallmark
sets2 <- read.gmt("TEPA_data/GOBP_CELL_MOTILITY.v2023.1.Mm.gmt")
sets3 <- read.gmt("TEPA_data/REACTOME_NEUTROPHIL_DEGRANULATION.v2022.1.Mm.gmt") # The only Reactome pathway we're interested in
sets4 <- read.gmt("TEPA_data/REACTOME_ENOS_ACTIVATION.v2023.1.Mm.gmt") # REACTOME_ENOS_ACTIVATION
sets5 <- read.gmt("TEPA_data/REACTOME_PENTOSE_PHOSPHATE_PATHWAY.v2023.1.Mm.gmt")
sets6 <- read.gmt("TEPA_data/WP_OXIDATIVE_STRESS_AND_REDOX_PATHWAY.v2023.1.Mm.gmt")
sets7 <- read.gmt("TEPA_data/WP_OXIDATIVE_STRESS_RESPONSE.v2023.1.Mm.gmt")
sets8 <- read.gmt("TEPA_data/WP_OXIDATIVE_DAMAGE_RESPONSE.v2023.1.Mm.gmt")
sets9 <- read.gmt("TEPA_data/REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE.v2023.1.Mm.gmt")
sets10 <- read.gmt("TEPA_data/REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE.v2023.1.Mm.gmt")

sets1$term <- as.character(sets1$term)
sets2$term <- as.character(sets2$term)
sets3$term <- as.character(sets3$term)
sets4$term <- as.character(sets4$term)
sets5$term <- as.character(sets5$term)
sets6$term <- as.character(sets6$term)
sets7$term <- as.character(sets7$term)
sets8$term <- as.character(sets8$term)
sets9$term <- as.character(sets9$term)
sets10$term <- as.character(sets10$term)

sets1 <- sets1 %>% split(x = .$gene, f = .$term)
sets2 <- sets2 %>% split(x = .$gene, f = .$term)
sets3 <- sets3 %>% split(x = .$gene, f = .$term)
sets4 <- sets4 %>% split(x = .$gene, f = .$term)
sets5 <- sets5 %>% split(x = .$gene, f = .$term)
sets6 <- sets6 %>% split(x = .$gene, f = .$term)
sets7 <- sets7 %>% split(x = .$gene, f = .$term)
sets8 <- sets8 %>% split(x = .$gene, f = .$term)
sets9 <- sets9 %>% split(x = .$gene, f = .$term)
sets10 <- sets10 %>% split(x = .$gene, f = .$term)
#sets3 <- sets3 %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets <- list()
fgsea_sets <- append(sets1, sets2)
fgsea_sets <- append(fgsea_sets, sets3)
fgsea_sets <- append(fgsea_sets, sets4)
fgsea_sets <- append(fgsea_sets, sets5)
fgsea_sets <- append(fgsea_sets, sets6)
fgsea_sets <- append(fgsea_sets, sets7)
fgsea_sets <- append(fgsea_sets, sets8)
fgsea_sets <- append(fgsea_sets, sets9)
fgsea_sets <- append(fgsea_sets, sets10)
fgsea_sets <- append(fgsea_sets, custom)
fgsea_sets <- append(fgsea_sets, copper_sign)

save = "N00_infTEPA_Enrichment"
Idents(seuset_nano) <- "Infiltration"
clusters = levels(Idents(seuset_nano))
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, out="pdf", input = "nanoInf_gCond")

### 7.3 - Treatment vs Control given Infiltration
save = "N00_TEPAInf_Enrichment_"
Idents(seuset_nano) <- "Condition"
clusters = levels(Idents(seuset_nano))
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, out="pdf", input = "nanoCond_gInf")

### 7.4 - Treatment vs Control 
save = "N00_TEPA_Enrichment_"
Idents(seuset_nano) <- "Condition"
clusters = levels(Idents(seuset_nano))
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, out="pdf", input = "nanoCond")

# TRY THE NETWORK


#### 8 - Deconvolution analysis using SingleR ####

#source("TEPA_code/DWLS-master/Deconvolution_functions.R")

# library(usethis) 
# usethis::edit_r_environ()

dataSC <- LoadH5Seurat("TEPA_results/S04_immuneDiff.h5Seurat")
Idents(dataSC) <- "scType"
# merge with tumor
dataBulk <- LoadH5Seurat("TEPA_results/N00_seusetNanoRed.h5Seurat", assay = "GeoMx")
labels = unique(Idents(dataSC))
scdata <- GetAssayData(dataSC[["integrated"]], slot = "scale.data")
spdata <- GetAssayData(seuset_nano, slot = "scale.data")


library(devtools)
library("SingleR")

pred <- SingleR(test=dataSC, ref=dataBulk, labels=labels, de.method="wilcox")
table(pred$labels)

### 9 - Expression of single genes and gene sets ####
# Search all isoforms of gene of interest
grep(pattern = "Cd8", 
     x = rownames(x = seuset_nano@assays$GeoMx@data), 
     value = TRUE, ignore.case = TRUE)

png("TEPA_plots/N00_Cd34_CompareCond.png", h = 2000, w = 3500, res = 200)
Idents(seuset_nano) <- "Condition"
DoHeatmap(object = subset(seuset_nano, downsample = 500), size = 6,
          features = c("Cd34", "Itgam", "Ptprc")) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + 
  theme(plot.margin = margin(2,2,1.5,1.2, "cm")) 
dev.off()

png("TEPA_plots/N00_Cd34_CompareCond.png", h = 2000, w = 3500, res = 200)
Idents(seuset_nano) <- "Condition"
DoHeatmap(object = subset(seuset_nano, downsample = 500), size = 6,
          features = c("Cd34", "Itgam", "Ptprc")) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + 
  theme(plot.margin = margin(2,2,1.5,1.2, "cm")) 
dev.off()

#png("TEPA_plots/S05_tumorMt1.png", h = 2000, w = 2000, res = 200)
pdf("TEPA_final_figures/N00_tumorMt1.pdf", h = 4, w = 4)
Idents(seuset_nano) <- "Infiltration"
levels(Idents(seuset_nano)) <- c("F", "T")
VlnPlot(seuset_nano, features =  c("Mt1"), assay = "GeoMx", slot="scale.data",
        cols = inf_col, ncol = 1, pt.size = 0) + 
  ylim(0,1.5) +
  #geom_signif(xmin = 1, xmax = 1.5, y_position = 2, annotations="***")+
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
  geom_signif(xmin = 1, xmax = 2, y_position = 5.25, annotations="***") +
  geom_boxplot(width=.1, outlier.size = 0.5)
dev.off()



immune_evasion <- c("Adora2a", "Arg1", "Arg2", "Cd274", "Cd276", "Cd47", "Cd80", "Cd83", "Cd86",
                    "Icosl", "Ido1", "Lgals3", "Pvr", "Tgfb1", "Tlr1", "Tlr2", "Tlr3", "Tlr4",
                    "Tlr5", "Tlr6", "Tlr7", "Tnfsf4", "Tnfsf9", "H2-D1", "H2-Ab1", "H2-Aa",
                    "H2-Eb1", "H2-Eb2", "H2-Ea", "H-2Dq", "H-2Lq")

my_data <- AverageExpression(
  seuset_nano,
  features = copper_genes,
  group.by = c("Infiltration","Condition"),
  slot = "scale.data")$GeoMx

# order of annotations/colors are defined here
ordered_meta_data <- str_split_fixed(colnames(my_data), '_', 2)
rownames(ordered_meta_data) <- colnames(my_data)   
colnames(ordered_meta_data) <- c("Infiltration", "Condition")

annotation_colors <- list("Infiltration" = inf_col,
                          "Condition" = cond_col)

ha = HeatmapAnnotation(df = as.data.frame(ordered_meta_data),
                       show_legend = TRUE,
                       show_annotation_name = TRUE,
                       col = annotation_colors)


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

save = "N00_complexHeat_Copper_sign"
#png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 15, w = 15)

Heatmap(
  my_data,
  col = col_fun,
  cluster_rows = TRUE,
  #row_km = ifelse(cluster,k,1),
  heatmap_legend_param=list(title="z-score"),
  row_names_gp = gpar(fontsize = 15, color = "white", lwd = 2),
  cluster_columns = FALSE,
  column_order = NULL,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_names_rot = 45,
  bottom_annotation = NULL,
  top_annotation = ha,
  use_raster = FALSE,
  heatmap_width = unit(25, "cm"), 
  heatmap_height = unit(25, "cm")
  #raster_by_magick = TRUE,
  #raster_quality = 5,
  #raster_resize_mat = TRUE
)
dev.off()


