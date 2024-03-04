#### Automatic clustering annotation with scType and DEA ##
## author: Antonietta Salerno
## date: 16/12/2022

library("Seurat")
library("writexl")
library(openxlsx)
library('limma')
library(dplyr)
library("MAST")
library(RColorBrewer)


setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadSeuratRds("TEPA_results/S02_immuneAnn.Rds")

#### 1 - Inter-cluster DEA: get marker genes ####

DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "celltypes"
clusters = levels(Idents(seuset_immune))
seuset_immune <- JoinLayers(seuset_immune)

# A - Find markers for every cluster compared to all remaining cells
immune.markers <- FindAllMarkers(seuset_immune, 
                                 only.pos = FALSE, 
                                 min.pct = 0.5, 
                                 min.diff.pct = 0.2,
                                 logfc.threshold = 0.5, 
                                 test.use="MAST",
                                 latent.vars="orig.ident")
immune.markers$p_val_adj = p.adjust(immune.markers$p_val, method='BH')
write.csv(immune.markers, "TEPA_results/S03_DEA_clusterMarkers.csv")

# Save results in different excel sheets 
clusters = levels(Idents(seuset_immune))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = immune.markers[immune.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/S03_DEA_clusterMarkers.xlsx", overwrite = TRUE)

# B - Add module score to the Seurat object for each cluster ###

seuset_immune <- createSets()

Idents(seuset_immune) <- "celltypes"
png("TEPA_plots/S03_immuneClustersAnnot.png", h = 3000, w = 4500, res = 200)
#pdf(qq("TEPA_final_figures/S03_immuneClustersAnnot.pdf"), h = 10, w = 14)
patchwork::wrap_plots(FeaturePlot(seuset_immune, ncol = 5, combine = TRUE,
                                  features = as.character(clusters), label = FALSE, repel = TRUE)) &
  theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off() 

# C - Identify most important markers per cluster ###

markers <- getTopMarkers(immune.markers, 5)

png("TEPA_plots/S03_immuneDotPlot.png", h = 2000, w = 2500, res = 300)
#pdf(qq("TEPA_final_figures/S03_immuneDotPlot.pdf"), h = 4, w = 11)
DotPlot(object = seuset_immune, features = unique(markers), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 4, cols = c("blue","red")) + RotatedAxis() + #scale_x_reverse() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=8))
dev.off()

png("TEPA_plots/S03_immuneMarkersHeatmap.png", h = 5000, w = 6000, res = 300)
top10 <- immune.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(object = subset(seuset_immune, downsample = 500), size = 6, 
          assay = "integrated", features = top10$gene) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

# D - Plot Volcano of each cluster vs all the others: 
save = "S03_immuneMarkers"
plotVolcano(clusters, res = immune.markers, type = "markers", immune = TRUE,
            log2FC = 0.5, save = save)

#### 2 - Intra-cluster DEA with annotated dataset - Treatment vs Control ####

seuset_immune$celltypes.tepa <- paste(Idents(seuset_immune), seuset_immune$condition, sep = "_")
seuset_immune$celltypes <- Idents(seuset_immune)
Idents(seuset_immune) <- "celltypes.tepa"
DefaultAssay(seuset_immune) <- "RNA"

save = "S03_immuneCond_"
sheets <- list()
for (cluster in unique(seuset_immune$celltypes)){
  try({
    ident1 <- paste0(cluster,"_Treatment")
    ident2 <- paste0(cluster,"_Control")
    condition.diffgenes <- FindMarkers(seuset_immune, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       #logfc.threshold = 0.2, 
                                       only.pos = FALSE, verbose = FALSE,
                                       latent.vars="orig.ident",
                                       #min.cells.feature = 1, min.cells.group = 1, 
                                       test.use="MAST")
    condition.diffgenes$p_val_adj = p.adjust(condition.diffgenes$p_val, method='BH')
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    write.csv(condition.diffgenes, file=paste0("TEPA_results/", save, "DEA_", gsub(" |/", "_", cluster),".csv"))
  })
}

# Needed for manual curation
openxlsx::write.xlsx(sheets, file = paste0("TEPA_results/", save, "DEA.xlsx"), rowNames=TRUE)

# Plot Volcano DEA by condition

Idents(seuset_immune) <- "celltypes"
clusters = unique(Idents(seuset_immune))

plotVolcano(clusters, log2FC = 0.25, pval = 0.05, save = save)

#### 3 - DEA Treatment vs Control bulk dataset ####

Idents(seuset_immune) <- "condition"
res <- FindMarkers(seuset_immune, ident.1 = "Treatment", ident.2 = "Control",
                                   only.pos = FALSE, verbose = FALSE,
                                   latent.vars="orig.ident",
                                   test.use="MAST")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/S03_immuneBulkDEA.csv"))

log2FC = 0.75
save = "S03_immuneBulk_"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     title = "DEA bulk TEPA vs Control",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 3,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")

SaveSeuratRds(seuset_immune,"TEPA_results/S03_immuneDiff.Rds")



