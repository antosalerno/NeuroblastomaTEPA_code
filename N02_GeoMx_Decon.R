#### Analysis GeoMx DSP data: deconvolution####
## author: Antonietta Salerno
## date: 17/06/2024
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

#### SpatialDecon ####

BiocManager::install("SpatialDecon")
library(SpatialDecon)

# Run SpatialDecon
results <- decon(
  dataset = dataset,
  segmentAnnotations = SampleAnnotationFile,
  cell_profile = aggregated_expression_matrix
)

# View the results
head(results)

#source("TEPA_code/DWLS-master/Deconvolution_functions.R")


#### SingleR ####
library(devtools)
library("SingleR")

dataSC <- LoadSeuratRds("TEPA_results/S08_seusetFull.SeuratRds")
Idents(dataSC) <- "scType"
# merge with tumor
dataBulk <- LoadSeuratRds("TEPA_results/N00_seusetNanoRed.Rds", assay = "GeoMx")
labels = unique(Idents(dataSC))
scdata <- GetAssayData(dataSC[["integrated"]], slot = "scale.data")
spdata <- GetAssayData(seuset_nano, slot = "scale.data")

pred <- SingleR(test=dataSC, ref=dataBulk, labels=labels, de.method="wilcox")
table(pred$labels)

library(usethis) 
usethis::edit_r_environ()
