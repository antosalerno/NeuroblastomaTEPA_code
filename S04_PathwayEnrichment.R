#### Pathway Enrichment Analysis ####
## author: Antonietta Salerno
## date: 19/12/2022

library("fgsea")
library(dplyr)
library(devtools)
library(msigdbr)
library("Seurat")
library(readr)
library(stringr)
library(clusterProfiler)
library(ggplot2)
library(ggcharts)
library("EnhancedVolcano")
library(GetoptLong)
library(magick)
library(tidyverse)
library(ComplexHeatmap)
# install.packages("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project/TEPA_data/org.Mm.eg.db_3.8.2.tar.gz") #download from source
library(org.Mm.eg.db)

#library(pheatmap)

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadSeuratRds("TEPA_results/S03_immuneDiff.Rds")
clusters = unique(seuset_immune@meta.data$celltypes)
immune.markers <- read.csv("TEPA_results/S03_DEA_clusterMarkers.csv")


#### 1 - Select the gene set collections of interest #### 
# HALLMARK_GLYCOLYSIS

# gene sets downloaded from source
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

#### 2 - Add custom gene set signatures ####
# Search all isoforms of gene of interest
grep(pattern = "Sigl", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

seuset_immune@assays$RNA@layers$scale.data <- scale(seuset_immune@assays$RNA@layers$data, scale = TRUE)
neutroCells = subset(seuset_immune, celltypes == "Neutrophils")

png("TEPA_plots/S04_neutrophilsPolarHeatmap.png", h = 4000, w = 6000, res = 300)
DoHeatmap(object = neutroCells, size = 6,slot="data",
          assay = "RNA", features = c(N1,N2), group.bar = T, group.colors = cond_col) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + 
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

# Plot Dotplot for the custom gene signatures

save = "S04_complexDot_Neutrophils_sign"
png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
#pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 15, w = 15)
sign_dotPlot(seuset_immune, c(N1,N2), cluster = FALSE, legend = FALSE) #check why it doesn't work
dev.off()

# Plot heatmap of average expression by celltype and condition

save = "S04_complexHeat_OnlyNeutrophils_sign"
png(paste0("TEPA_plots/", save, ".png"), h = 4000, w = 2500, res = 300)
#pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 10, w = 5)
neutroCells <- subset(seuset_immune, celltypes == "Neutrophils")
sign_avgHeatMap(neutroCells, c(N1,N2), cluster = FALSE, legend=FALSE, w=5, h=25)
dev.off()

# Plot signature on the UMAP

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(N1), name=make.names("N1"))
names(seuset_immune@meta.data)[grep(make.names("N1"), names(seuset_immune@meta.data))] <- "N1"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(N2), name=make.names("N2"))
names(seuset_immune@meta.data)[grep(make.names("N2"), names(seuset_immune@meta.data))] <- "N2"

png("TEPA_plots/S04_N1_N2_FeaturePlot.png", h = 2000, w = 2500, res = 250)
#pdf(qq("TEPA_final_figures/S04_N1_N2_FeaturePlot.pdf"), h = 8, w = 10)
FeaturePlot(seuset_immune, ncol = 2, pt.size = 0.5, split.by = "condition",
            features = c("N1", "N2"), label = F,
            repel = TRUE) & theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

#### 3 - Run the custom gsea function ####
clusters = unique(levels(seuset_immune$celltypes))
save = "S04_immuneEnrichment_"
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, minSize = 6, adj = TRUE, out = "pdf")

#gseaPlotRes(clusters)

#### 4 - Plot all clusters' results in a network ####

save = "S04_immuneJointNet"
gseaByCellType <- gseaJointNet(clusters, save = save)

# If you want to filter out some gene set manually!
df_gseaCT <- read.csv(paste0("TEPA_results/", save, "SHORT3.csv"), sep = ";")
gseaByCellType@compareClusterResult <- df_gseaCT

cnet <- cnetplot(gseaByCellType, showCategory=6, 
                 color_category = cellt_col[1:13],
                 color.params = list(category = cellt_col[1:13]),
                 cex.params = list(gene_label = 2, category_label = 3)) 
file=paste0("TEPA_plots/",save,".pdf")
ggsave(cnet, file = file, width = 55, height = 45, units = "cm")


#### 5 - Clustered diverging bar plot all cell types

save = "S04_immuneJointBarplot3"
gseaByType(clusters, save = save)

fgseaResByType = read.csv(paste0("TEPA_results/", save, "SHORT.csv"), sep = ";")
b <- barPlotGSEA(fgseaResByType, byType = TRUE)
ggsave(b, file=paste0("TEPA_plots/S04_barplotCellTypesEnriched3.pdf"),
       width = 40, height = 20, units = "cm", limitsize = F, dpi = 500)

SaveSeuratRds(seuset_immune, "TEPA_results/S04_immuneDiff.Rds")

