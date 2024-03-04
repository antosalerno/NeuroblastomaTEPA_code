###### Dimensionality reduction and clustering of immune cells with higher resolution##
## author: Antonietta Salerno
## date: 13/01/2022

library("Seurat")
library("ggplot2")
library("writexl")
library(openxlsx)
library(HGNChelper)
library(dplyr)
library(RColorBrewer)
library(patchwork)

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")

#### 1 - Dimensionality reduction ####

immune.combined <- LoadSeuratRds("TEPA_results/S01_immuneInt.Rds")

# Scaling
seuset_immune <- ScaleData(immune.combined, verbose = FALSE)
seuset_immune <- RunPCA(seuset_immune, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_immune@reductions$pca@stdev / sum(seuset_immune@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 70.26% of variability

seuset_immune <- FindNeighbors(object = seuset_immune, graph.name = "clust", dims = 1:60, reduction = 'pca')

### 2 - Clustering ####

seuset_immune <- FindClusters(object = seuset_immune, graph.name = "clust", resolution = 0.3) 
seuset_immune <- RunUMAP(seuset_immune, dims = 1:60, reduction = "pca", verbose = FALSE)

png("TEPA_plots/S02_umapExplore.png", w = 6000, h = 2000, res = 300)
#pdf(qq("TEPA_final_figures/S02_umapExplore.pdf"), w = 15, h = 8)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("orig.ident", "condition", "seurat_clusters"), label = FALSE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("TEPA_plots/S02_umapCounts.png", w = 4000, h = 2000, res = 300)
FeaturePlot(seuset_immune, features = c("nCount_RNA", "nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

png("TEPA_plots/S02_umapClust.png", w = 4000, h = 2000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("seurat_clusters"), split.by= "condition",label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Sub-clustering ##

Idents(seuset_immune) <- "clust_res.0.3"
seuset_immune <- FindSubCluster(seuset_immune, "3", "clust", subcluster.name = "seurat_clusters",  resolution = 0.2, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"

seuset_immune <- FindSubCluster(seuset_immune, "3_0", "clust", subcluster.name = "seurat_clusters",  resolution = 0.3, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"

seuset_immune <- FindSubCluster(seuset_immune, "2", "clust", subcluster.name = "seurat_clusters",  resolution = 0.35, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"

seuset_immune <- FindSubCluster(seuset_immune, "2_2", "clust", subcluster.name = "seurat_clusters",  resolution = 0.2, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"


### 4 - Annotate subclusters
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Immune system"

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = seuset_immune[["integrated"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

####  4 - Get scType scores by cluster ####

cL_results = do.call("rbind", 
                     lapply(unique(seuset_immune@meta.data$seurat_clusters), 
                            function(cl){
                              es.max.cl = sort(rowSums(es.max[,rownames(seuset_immune@meta.data[seuset_immune@meta.data$seurat_clusters==cl, ])]),decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuset_immune@meta.data$seurat_clusters==cl)), 10)
                            }))
celltypes_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
celltypes_scores = celltypes_scores[order(celltypes_scores$cluster),]

# Manual curation of celltypes annotation
celltypes_scores[3, "type"] <- "Cd4+ Naive T cells"
celltypes_scores[4, "type"] <- "Cd4+ Memory T cells" #"Cd4+ Naive T cells"
celltypes_scores[5, "type"] <- "Cd8+ Naive-Memory T cells" #"Cd4+ Memory T cells"
celltypes_scores[6, "type"] <- "Cd8+ NkT-like cells"
celltypes_scores[7, "type"] <- "DN Regulatory T cells" #"Cd8+ Naive T cells"
celltypes_scores[8, "type"] <- "Cd4+ Naive T cells" #"Cd8+ NkT-like cells"
celltypes_scores[10, "type"] <- "Macrophages" # or macrophages
celltypes_scores[11, "type"] <- "Macrophages" # or macrophages
celltypes_scores[12, "type"] <- "HLA-expressing cells" #"Dendritic cells" # or DCs
celltypes_scores[13, "type"] <- "Dendritic cells" #"B cells"
celltypes_scores[14,"type"] <- "Natural killer cells"
celltypes_scores[15, "type"] <- "B cells"
celltypes_scores[17,"type"] <- "Gamma-delta T cells"

seuset_immune@meta.data$celltypes = ""
for(j in unique(celltypes_scores$cluster)){
  cl_type = celltypes_scores[celltypes_scores$cluster==j,]; 
  seuset_immune@meta.data$celltypes[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

Idents(seuset_immune) <- seuset_immune@meta.data$celltypes


# We decide to remove the cluster HLA-expressing cells since it is probably ambient RNA -> mixture of B cells and macrophage markers
seuset_immune <- seuset_immune[,!seuset_immune$celltypes == "HLA-expressing cells"] # 12 cells
# We are not sure how to define these progenitors
seuset_immune <- seuset_immune[,!seuset_immune$celltypes == "Progenitor cells"] # 45 cells

levels(Idents(seuset_immune)) # now 14 clusters rather than 15

seuset_immune@meta.data$celltypes <- factor(seuset_immune@meta.data$celltypes,
                                  levels=c("Cd4+ Naive T cells","Cd4+ Memory T cells",  
                                           "Cd8+ Naive-Memory T cells","Cd8+ NkT-like cells" , 
                                           "Gamma-delta T cells","DN Regulatory T cells" , 
                                           "Natural killer cells", 
                                           "B cells" , "Dendritic cells", "Macrophages", 
                                           "Basophils", "Eosinophils", "Neutrophils"  
                                  ))
Idents(seuset_immune) <- "celltypes"
options(ggrepel.max.overlaps = Inf)
png("TEPA_plots/S02_umapAnn.png", w = 3000, h = 3000, res = 300)
p <- DimPlot(object = seuset_immune, pt.size = 0.5, reduction = 'umap', ncol = 1,
             group.by = "celltypes", label = FALSE,cols = cellt_col) +
  ggtitle("Cell types in NB Control and Treated samples") +
  theme(plot.title = element_text(hjust = 0.5)) 
LabelClusters(p, id = "celltypes", size = 5, repel = T,  box.padding = 1)
dev.off()

png("TEPA_plots/S02_umapCondAnn.png", w = 4000, h = 2000, res = 300)
#pdf(qq("TEPA_final_figures/S02_umapCondAnn.pdf"), w = 12, h = 6)
pt1 <- DimPlot(object = seuset_immune, pt.size = 0.08, reduction = 'umap', ncol = 2,
        group.by = "celltypes", split.by= "condition", label = FALSE, cols = cellt_col) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) 
LabelClusters(pt1, id = "celltypes", size = 2.8, repel = T)
dev.off()

### Bar Plot with proportions
pt <- table(Idents(seuset_immune), seuset_immune$condition)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1,levels=levels(seuset_immune$celltypes)) # Look the difference in the factors order

png(paste0("TEPA_plots/S02_condAnnClusterFreq.png"), w=2500,h=2500, res=300)
#pdf(qq("TEPA_final_figures/S02_condAnnClusterFreq.pdf"), w = 6, h = 8)
ggplot(pt, aes(x = Var2, y = Freq, fill = factor(Var1))) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity")+
  xlab("Condition") +
  ylab("Cell type") +
  #ggtitle("Cell types in TEPA vs Control") +
  scale_fill_manual(values = cellt_col) +
  theme(legend.title = element_blank()) & NoLegend()
dev.off()

###  Bar Plot with percentage

pt <- table(Idents(seuset_immune), seuset_immune$condition)
data_percentage<- apply(pt, 2, function(x){as.numeric(x)*100/sum(x,na.rm=T)})

# Make a stacked barplot--> it will be in %!
png(paste0("TEPA_plots/S02_condAnnClusterPerc.png"), w=2500,h=2500, res=300)
#pdf(paste0("TEPA_final_figures/S02_condAnnClusterPerc.pdf"), w=2500,h=2500)
par(mar = c(5.1, 5.1, 4.1, 12))
barplot(data_percentage, col=cellt_col, border="white", 
        ylab="Percentage of cells per cell type",
        main = "Cell types in TEPA vs Control",
        legend = rownames(pt),
        args.legend = list(x = "topright",inset = c(-0.45, 0)))
dev.off()

SaveSeuratRds(seuset_immune, "TEPA_results/S02_immuneAnn.SeuratRds")


