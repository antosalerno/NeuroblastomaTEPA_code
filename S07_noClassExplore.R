###### Explore the cells not classified as tumor nor immune ##
## author: Antonietta Salerno
## date: 04/03/2024

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


seuset <- LoadSeuratRds("TEPA_results/S00_seusetClass.Rds")
