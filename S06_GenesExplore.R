#### Validate clustering annotation ####
## author: Antonietta Salerno
## date: 20/12/2022

library("Seurat")
library("ggplot2")
library(RColorBrewer)
library("ggpubr")

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadSeuratRds("TEPA_results/S03_immuneDiff.Rds")
seuset_full<- LoadSeuratRds("TEPA_results/S08_seusetFull.Rds")
immune.markers <- read.csv("TEPA_results/02_DEA_clusterMarkers.csv")

DefaultAssay(seuset_immune) <- "RNA"

# Search all isoforms of gene of interest
grep(pattern = "H", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

#  Plot density plot with key markers for important populations  

png("TEPA_plots/03_densityPlatelets.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <-"RNA"
plot_density(seuset_immune, features = c("Itga2b", "Itgb3"), joint = TRUE)
dev.off()

png("TEPA_plots/03_densityErythro.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = "Gypa")
dev.off()

png("TEPA_plots/03_densityTreg.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Il2ra","Foxp3"), joint = TRUE)
dev.off() # Treg higher in tumor and in inflammation: cd25 and foxp3

# Classical monocytes and Nks: Ccr2 important for neutrophils mobilitisation and recruitment to metastatic sites

png("TEPA_plots/03_densityEosino.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Il5ra", "Ccr3", "Adgre1","Ccl6"), joint = TRUE)
dev.off()

png("TEPA_plots/03_densityMemory.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Il7r", "Ctla4"), joint = TRUE)
dev.off()

png("TEPA_plots/03_Mt2_split.png", h = 2000, w = 3500, res = 200)
plot_density(seuset_immune, "Mt2") + 
  facet_grid(.~seuset_immune$condition)
#plot_density(seuset_immune, features = c("Mt1", "Mt2"), joint = TRUE)
dev.off()

png("TEPA_plots/03_Mt1Mt2_violin.png", h = 2000, w = 3500, res = 200)
Idents(seuset_immune) <- "condition"
VlnPlot(seuset_immune, features = c("Mt2", "Mt1"), ncol = 2,pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_Ifngr_violin.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
VlnPlot(seuset_immune, features = "Ifngr1", split.by="condition", ncol = 1,pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_Ceacam.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
VlnPlot(seuset_immune, features = c("Ceacam1","Ceacam10"), 
        split.by="condition", ncol = 2, pt.size = 0.000005) + geom_boxplot()
dev.off()

png("TEPA_plots/03_NeutrophilsActive.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Fut4", "Itgb2", "Itgam"), 
        split.by="condition", ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_NeutroMacro.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Ly6g", "Ly6g5b", "Ly6g6c"), ncol = 3, pt.size = 0.000005)
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


#### See difference in proportions of lymphoid and myeloid cells tumor vs immune cells ####

png("TEPA_plots/S06_LymphoMyeloCompare.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full) <- "class"
VlnPlot(seuset_full, features = c("Cd3e", "Itgam"), split.by = "condition",
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

#### Check Control T-cells ####

grep(pattern = "Th1", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)


CD8_exh <- c("Pdcd1", "Havcr2", "Lag3", "Ctla4", "Cd244a", "Cd160")
#PD-1, TIM-3, LAG3, CTLA-4, 2B4, Cd244, TIGIT
CD4_reg <- c("Foxp3", "Il2ra", "Il10")
# FOXP3, CD25, Il10


supp <- c(CD8_exh, CD4_reg)
save = "S06_complexHeat_Exh"
png(paste0("TEPA_plots/",save,".png"), h = 5000, w = 6000, res = 400)
sign_avgHeatMap(seuset_full, CD8_exh, immune = FALSE,
                cluster = TRUE, k = 3, legend = TRUE) #check why it doesn't work
dev.off()

immune_sup_cancer <- c("Ido1", "Cd274") #Pdl1


png("TEPA_plots/S06_Ido1_VlnPlot.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full) <- "class"
VlnPlot(seuset_full, features = immune_sup_cancer, split.by = "condition",
        ncol = 2, pt.size = 0.000005)
dev.off()



