#### Analysis GeoMx DSP data: Enrichment ####
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

#### 7 - Pathway enrichment Analysis ####
seuset_nano <- LoadSeuratRds("TEPA_results/N00_seusetNanoRed.Rds")

### 7.1 - Cluster markers
### 7.2 - Infiltration vs Non-Infiltration given Treatment

library(dplyr)

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
pdf("TEPA_final_figures/N00_tumorMts.pdf", h = 4, w = 4)
Idents(seuset_nano) <- "Infiltration"
levels(Idents(seuset_nano)) <- c("F", "T")
VlnPlot(seuset_nano, features =  c("Mt1", "Mt2"), assay = "GeoMx", slot="scale.data",
        cols = inf_col, ncol = 3, pt.size = 0) + 
  ylim(0,1.5) +
  #geom_signif(xmin = 1, xmax = 1.5, y_position = 2, annotations="***")+
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
  layer = "scale.data")$GeoMx

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

save = "N00_complexHeat_Copper_sign2"
#png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 15)

Heatmap(
  as.matrix(my_data),
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


