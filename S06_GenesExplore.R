#### Validate clustering annotation ####
## author: Antonietta Salerno
## date: 20/12/2022

library("Seurat")
library("ggplot2")
library(RColorBrewer)
library("ggpubr")

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
#seuset_immune <- LoadSeuratRds("TEPA_results/S03_immuneDiff.SeuratRds")
#seuset_full<- LoadSeuratRds("TEPA_results/S08_seusetFull.SeuratRds")
immune.markers <- read.csv("TEPA_results/S03_DEA_clusterMarkers.csv")
immune.markers <- read.csv("TEPA_results/S03_immuneCond_DEA.xlsx")


DefaultAssay(seuset_immune) <- "RNA"

# Search all isoforms of gene of interest
grep(pattern = "H", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

png("TEPA_plots/03_Mt2_split.png", h = 2000, w = 3500, res = 200)
plot_density(seuset_immune, "Mt2") + 
  facet_grid(.~seuset_immune$condition)
#plot_density(seuset_immune, features = c("Mt1", "Mt2"), joint = TRUE)
dev.off()

png("TEPA_plots/03_Mt1Mt2_violin.png", h = 2000, w = 3500, res = 200)
Idents(seuset_immune) <- "condition"
VlnPlot(seuset_immune, features = c("Mt2", "Mt1"), ncol = 2,pt.size = 0.000005)
dev.off()

png("TEPA_plots/S03_Prnp_ImmuneCond.png", h = 2000, w = 3500, res = 400)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
VlnPlot(seuset_immune, features = c("Prnp"), 
        split.by = "condition", ncol = 1, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S03_Mt1_ImmuneCond.png", h = 2000, w = 3500, res = 400)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
fig <- VlnPlot(seuset_immune, features = c("Mt1"), 
               split.by = "condition", ncol = 1, pt.size = 0.000005)
fig +  ylim(0, 5) +
  geom_signif(xmin = 0.75, xmax = 1.25, y_position = 3.75, annotations="***") +
  geom_signif(xmin = 4.75, xmax = 5.25, y_position = 4.5, annotations="***")
dev.off()


#### See difference in proportions of lymphoid and myeloid cells tumor vs immune cells ####

png("TEPA_plots/S06_LymphoMyeloCompare.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full) <- "class"
VlnPlot(seuset_full, features = c("Cd3e", "Itgam"), split.by = "condition",
            ncol = 2, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S03_lympho.png", h = 4000, w = 6500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Cd4", "Cd8a", "Cd8b1", "Cd3d", "Cd3e", "Cd3g"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S06_LymphoMyelo.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Cd3e", "Itgam"),
            ncol = 2, pt.size = 0.000005)
dev.off()

### Check pathway of immune evasion ####

immune_evasion <- c("Adora2a", "Arg1", "Arg2", "Cd274", "Cd276", "Cd47", "Cd80", "Cd83", "Cd86",
                    "Icosl", "Ido1", "Lgals3", "Pvr", "Tgfb1", "Tlr1", "Tlr2", "Tlr3", "Tlr4",
                    "Tlr5", "Tlr6", "Tlr7", "Tnfsf4", "Tnfsf9", "H2.D1", "H2.Ab1", "H2.Aa",
                    "H2.Eb1", "H2.Eb2", "H2.Ea", "H.2Dq", "H.2Lq")


save = "S06_complexHeat_ImmuneEvasion"
png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
#pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 15, w = 15)
sign_avgHeatMap(seuset_full, immune_evasion, immune = FALSE,
                cluster = FALSE, k = 2, legend = TRUE) #check why it doesn't work
dev.off()

### Check cytokines assay ####

cytokines <- c("Csf2","Csf3", "Il10", "Lif", "Il1b", "Il2","Csf1", "Cxcl10",
               "Il4", "Il5","Il6", "Ifnar1", "Ifnar2", "Il3ra", "Il22", "Il13", "Il27",
               "Il23a", "Ifng", "Il12a", "Il12b", "Cxcl1", "Ccl5", "Tnf",
               "Ccl3", "Ccl7", "Ccl2", "Il17a", "Il15", "Cxcl2", "Il1a", 
                "Ccl11", "Il18", "Ccl4", "Ifnb1", "Ifnlr1", "Il9r")

save = "S06_complexHeat_Cytokines"
#png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
pdf(qq(paste0("TEPA_final_figures/",save,".pdf")), h = 15, w = 15)
sign_avgHeatMap(seuset_full, cytokines, immune = FALSE,
                cluster = TRUE, k = 3, legend = TRUE) #check why it doesn't work
dev.off()

#### Immune activation and migration ####

immune_migr <- c("Trem1", "Adam8", "Clec4e", "Cd14", "Ccr7", "Cxcr4", "Fpr1")

immune_extra <-c("Itgam", "Itgb2", "Cd177", "Itgal", "Cd44", "Itga2", "Itgb1", "Myo1f")

interferon <-c("Ifit1","Ifit2", "Ifit3", "Isg15", "Isg20", "Ifitm1", "Ifitm2", "Ifitm3", "Ifitm6")


save = "S06_complexHeat_Migration"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 15)
sign_avgHeatMap(seuset_full, immune_migr, immune = FALSE,
                cluster = TRUE, k = 1, legend = TRUE) 
dev.off()

save = "S06_complexHeat_Extravasation_scaled"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 15)
sign_avgHeatMap(seuset_full, immune_extra, immune = FALSE,
                cluster = TRUE, k = 1, legend = TRUE) 
dev.off()

save = "S06_complexHeat_ExtraMigr_FC_Neutro"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 5)
sign_avgLogFCHeatMap(seuset_full, c(immune_extra, immune_migr), immune = TRUE, 
                     celltype = "Neutrophils",
                cluster = TRUE, k = 1, legend = TRUE) 
dev.off()

save = "S06_complexHeat_Interferon_FC_Neutro"
pdf(paste0("TEPA_final_figures/",save,".pdf"), h = 15, w = 5)
sign_avgLogFCHeatMap(seuset_full, interferon, immune = FALSE, celltype = "Neutrophils",
                cluster = TRUE, k = 1, legend = TRUE) 
dev.off()




