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
library("fgsea")
library("stringr")

### 1 - Load data ####
setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
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
SaveSeuratRds(seuset_nano, "TEPA_results/N00_seusetNano.Rds")

#head(seuset_nano@misc$QCMetrics$QCFlags) 

png("TEPA_plots/N00_countsROIs.png", h = 3000, w = 2500, res = 300)
VlnPlot(seuset_nano, features = "nCount_GeoMx", split.by = "Infiltration",
        pt.size = 0.1)
dev.off() # the number of genes is instead the same for all ROIs

png("TEPA_plots/N00_nucleiROIs.png", h = 3000, w = 2500, res = 300)
VlnPlot(seuset_nano, features = "AOINucleiCount", split.by = "Group",
        pt.size = 0.1)
dev.off() 

#### 5 - Dimensionality reduction ####
seuset_nano <- LoadSeuratRds("TEPA_results/N00_seusetNano.Rds")

seuset_nano <- FindVariableFeatures(seuset_nano)
seuset_nano <- ScaleData(seuset_nano)
seuset_nano <- RunPCA(seuset_nano, assay = "GeoMx", verbose = FALSE, approx=FALSE)

# Determine percent of variation associated with each PC
pct <- seuset_nano@reductions$pca@stdev / sum(seuset_nano@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=50) # Select 50 PCs to retain 99% of variability

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
seuset_immune <- LoadSeuratRds("TEPA_results/S04_immuneDiff.Rds")


# 5.2 Infiltration vs Non-Infiltration given Treatment

Idents(seuset_nano) <- "Infiltration"
seuset_nanoTEPA <- subset(seuset_nano, Condition == "Treatment")

seuset_nanoTEPA@assays$GeoMx@layers$scale.data <- scale(seuset_nanoTEPA@assays$GeoMx@layers$counts)

seuset_nanoTEPA <- FindVariableFeatures(seuset_nanoTEPA, selection.method = "vst")

res <- FindMarkers(seuset_nanoTEPA, ident.1 = "T", ident.2 = "F", slot="counts",
                   only.pos = FALSE, verbose = FALSE, assay= "GeoMx", 
                   test.use="negbinom")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoInf_gCond_DEA.csv"))


SaveSeuratRds(seuset_nano, "TEPA_results/N00_seusetNanoRed.Rds")

