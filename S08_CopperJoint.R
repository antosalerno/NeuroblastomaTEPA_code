#### Joint analysis to understand copper redistribution in the immune cells after chelation ####
## author: Antonietta Salerno
## date: 13/04/2023

library("Seurat")
library("ggplot2")
library(RColorBrewer)
library(magick)
library(GetoptLong)


setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadSeuratRds("TEPA_results/S04_immuneDiff.Rds")
seuset_tumor <- LoadSeuratRds("TEPA_results/S05_seusetTumorClu.Rds")

seuset_full <- merge(seuset_immune, y = seuset_tumor, 
                     add.cell.ids = c("immune", "tumor"), 
                     project = "singleCell")

#seuset_full@assays$RNA@layers$scale.data <- scale(seuset_full@assays$RNA@layers$data, scale = TRUE)

seuset_full$celltypes <- ifelse(test = is.na(seuset_full$celltypes), yes = "Tumor", no = seuset_full$celltypes)

seuset_full$condition <- factor(seuset_full$condition,
                                levels=c("Control", "Treatment"))

seuset_full@meta.data$celltypes <- factor(seuset_full@meta.data$celltypes,
                                            levels=c("Cd4+ Naive T cells","Cd4+ Memory T cells",  
                                                     "Cd8+ Naive-Memory T cells","Cd8+ NkT-like cells" , 
                                                     "Gamma-delta T cells","DN Regulatory T cells" , 
                                                     "Natural killer cells", 
                                                     "B cells" , "Dendritic cells", "Macrophages", 
                                                     "Basophils", "Eosinophils", "Neutrophils", "Tumor"  
                                            ))

SaveSeuratRds(seuset_full, "TEPA_results/S08_seusetFull.Rds")


#### Create a signature of copper-related genes ####

copper_genes <- c("Sod1", "Sod2", "Gls","Sp1", "Atox1",  "Mtf2","Pdha1", "Pdhb", "Lias", "Dld", "Dlat", 
                  "Mt1", "Mt2", "Slc31a1", "Atp7a","Steap4",
                  "Atp7b", "Steap3", 
                  "Sco1", "Cox11", "Commd1", "Mtf1","Fdx1")


# Search all isoforms of gene of interest
grep(pattern = "Sod1", 
     x = rownames(x = seuset_full@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

#### Check differentially expressed genes in different cell types ####

seuset_full <- LoadSeuratRds("TEPA_results/S08_seusetFull.Rds")

clusters = levels(seuset_immune$celltypes)
df <- data.frame()
for (cluster in clusters){
  file=paste0("TEPA_results/S03_immuneCond_DEA_",gsub(" |/", "_", cluster),".csv")
  res <- read.csv(file, sep=",")
  res$celltype <- cluster
  df <- rbind(df, res)
}
file=paste0("TEPA_results/S05_tumorBulkDEA.csv")
res <- read.csv(file, sep=",")
res$celltype <- "Tumor"
df <- rbind(df, res)

df_copper <- df[df$X %in% copper_genes,] 
df_copper <- df_copper %>% filter(p_val < 0.05)

View(df_copper)

png("TEPA_plots/S08_violinCopperGenesSign.png", h = 6000, w = 2000)
Idents(seuset_full) <- "celltypes"
VlnPlot(seuset_full, features = unique(df_copper$X), 
        assay = "RNA", split.by = "condition",
        #size.x.use = 10, size.y.use = 10, 
        ncol = 1, pt.size = 0.000005) 
dev.off()


#### Visualise expression of copper genes ####

### A - Dotplot to understand what most of the cells in the immune subset express ###

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

save = "S08_complexDot_Copper_sign"
#png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 15, w = 15)
sign_dotPlot(seuset_immune, copper_genes, immune=TRUE, cluster = FALSE,k=2, legend = FALSE) #check why it doesn't work
dev.off()

#### B - Create an heatmap to understand differences between immune cells and tumor cells as well as control and treatment ####

save = "S08_complexHeat_Copper_sign"
#png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 15, w = 15)
sign_avgHeatMap(seuset_full, copper_genes, immune = FALSE,
                cluster = FALSE, k = 1, legend = FALSE) #check why it doesn't work
dev.off()


# The expression values for each gene are scaled / standardized by subtracting the genes mean expression and dividing by its standard deviation. 
# A value of -1 would imply it's one standard deviation below the mean expression for that gene.

### SINGLE-CELL ###

# order of annotations/colors are defined here
ordered_meta_data <- seuset_full@meta.data[order(seuset_full@meta.data$celltypes), ]
ordered_meta_data <- select(ordered_meta_data, c('condition','celltypes'))

annotation_colors <- list("celltypes" = cellt_col,
                          "condition" = cond_col)

ha = HeatmapAnnotation(df = as.data.frame(ordered_meta_data),
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# Expression data
my_data <- seuset_full[copper_genes,]@assays$RNA@scale.data
my_data <- my_data[, rownames(ordered_meta_data)]

# Heatmap
png("TEPA_plots/S08_fullGenesCopper_complexHeat_SingleCell.png", h = 2000, w = 2500, res = 300)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
h2 <- Heatmap(
  my_data,
  col = col_fun,
  cluster_rows = TRUE,
  heatmap_legend_param=list(title="z-score"),
  cluster_columns = FALSE,
  column_order = NULL,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  use_raster = FALSE,
  column_names_rot = 45,
  #raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha
)
dev.off()








