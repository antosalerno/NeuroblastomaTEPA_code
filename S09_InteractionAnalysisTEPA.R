#### Deciphering cell–cell interactions and communication from gene expression ####
## author: Antonietta Salerno
## date: 02/01/2024

library("Seurat")
library("ggplot2")
library(RColorBrewer)
library(magick)
library(GetoptLong)
library("NMF")
library("circlize")
library(ComplexHeatmap)
#devtools::install_github("jinworks/CellChat")
library("CellChat")
library(patchwork)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
seuset_full <- LoadSeuratRds("TEPA_results/S08_seusetFull.Rds")

# 1. DATA PREPARATION ####
seuset_tepa <- seuset_full[,seuset_full$condition == "Treatment"]

data.input <- GetAssayData(seuset_tepa, assay = "RNA", layer = "data") # normalized data matrix
labels <- seuset_tepa$celltypes
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchatTEPA <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchatTEPA <- addMeta(cellchatTEPA, meta = meta)
cellchatTEPA <- setIdent(cellchatTEPA, ident.use = "group") # set "labels" as default cell identity
levels(cellchatTEPA@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchatTEPA@idents)) # number of cells in each cell group

# 2. SET THE LIGAND-RECEPTOR INTERACTION DATABASE for cell communication analysis ####

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchatTEPA@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchatTEPA <- subsetData(cellchatTEPA) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchatTEPA <- identifyOverExpressedGenes(cellchatTEPA)
cellchatTEPA <- identifyOverExpressedInteractions(cellchatTEPA)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# 3. INFERENCE OF CELL-CELL COMMUNICATION METHOD ####

# Compute the communication probability and infer cellular communication network
cellchatTEPA <- computeCommunProb(cellchatTEPA, type = "triMean", population.size = TRUE)
#> triMean is used for calculating the average gene expression per cell group. 
#> ‘trimean’ approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. >> most important interactions very stringent
#> To use 10% truncated mean, USER can set type = "truncatedMean" and trim = 0.1. 
#> To determine a proper value of trim, function computeAveExpr can check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1)

cellchatTEPA <- filterCommunication(cellchatTEPA, min.cells = 10)

df.net <- subsetCommunication(cellchatTEPA)

df.net_short <- subsetCommunication(cellchatTEPA, sources.use = c("Neutrophils"), targets.use = c("Gamma-delta T cells","Cd4+ Naive T cells", "DN Regulatory T cells",
                                                                                              "Cd8+ Naive-Memory T cells", "Cd8+ NkT-like cells", "Natural killer cells",
                                                                                              "Macrophages", "B cells", "Basophils", "Eosinophils", "Dendritic cells", "Tumor"))

#> df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

cellchatTEPA <- computeCommunProbPathway(cellchatTEPA)

cellchatTEPA <- aggregateNet(cellchatTEPA)

groupSize <- as.numeric(table(cellchatTEPA@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("TEPA_project/TEPA-figures/S09_TEPAcellComm1.pdf", h = 10, w = 10)
netVisual_circle(cellchatTEPA@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatTEPA@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("TEPA_project/TEPA-figures/S09_TEPAcellComm2.pdf", h = 10, w = 15)
mat <- cellchatTEPA@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# 4. VISUALISATION ####

# A- Hierarchy plot
pathways.show <- c("TNF") 
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4)# a numeric vector. 
par(mfrow = c(1,1), xpd=TRUE)
pdf("TEPA_final_figures/S09_TEPAcellComm_TNF.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_aggregate(cellchatTEPA, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# B- Chord diagram -> does it take into account population proportions?
par(mfrow=c(1,1))
pathways.show <- c("CXCL") 
pdf("TEPA_final_figures/S09_TEPAcellComm_CXCL_chord.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_aggregate(cellchatTEPA, signaling = pathways.show, layout = "chord")
dev.off()

par(mfrow=c(1,1))
pathways.show <- c("TNF") 
pdf("TEPA_final_figures/S09_TEPAcellComm_TNF_chord.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_aggregate(cellchatTEPA, signaling = pathways.show, layout = "chord")
dev.off()

# C - Heatmap
par(mfrow=c(1,1))
pathways.show <- c("CXCL") 
pdf("TEPA_final_figures/S09_TEPAcellComm_TNF_heatmap.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_heatmap(cellchatTEPA, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

# see contribution of each LR pair
pathways.show <- c("TNF")  # Tnfrsf1b/a
netAnalysis_contribution(cellchatTEPA, signaling = pathways.show) # you can include only one LR pair in the analysis

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchatTEPA@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchatTEPA@idents)
vertex.receiver = seq(1,4)
pdf("TEPA_final_figures/S09_TEPAcellComm_All_chord.pdf", h = 10, w = 10)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_aggregate(cellchatTEPA, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "chord")
}
dev.off()

pdf("TEPA_final_figures/S09_TEPAcellComm_All_circle.pdf", h = 12, w = 10)
for (i in 1:length(pathways.show.all)) {
  netVisual_aggregate(cellchatTEPA, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver,
                      signaling.name = pathways.show.all[i])
}
dev.off()

# Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
for (i in 1:length(pathways.show.all)) {
  p<-netAnalysis_contribution(cellchatTEPA, signaling = pathways.show.all[i])
  file = paste0("TEPA_final_figures/S09_TEPAcellComm_contribution_", pathways.show.all[i])
  ggsave(p, filename = file, width = 4, height = 2, units = 'in', dpi = 300, device = "pdf")
}

# D- Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
levels(cellchatTEPA@meta$group)
pdf("TEPA_project/TEPA-figures/S09_TEPAcellComm_neutrophils_vs_all_BUBBLE.pdf", h = 6, w = 5)
netVisual_bubble(cellchatTEPA, sources.use = "Neutrophils", targets.use = c(1:14), remove.isolate = FALSE)
dev.off()

levels(cellchatTEPA@meta$group)
pdf("TEPA_project/TEPA-figures/S09_TEPAcellComm_all_vs_neutrophils_BUBBLE.pdf", h = 6, w = 5)
netVisual_bubble(cellchatTEPA, sources.use = c(1:14), targets.use = "Neutrophils", remove.isolate = FALSE)
dev.off()

# 5. Systems analysis of cell-cell communication network ####

# Compute the network centrality scores
cellchatTEPA <- netAnalysis_computeCentrality(cellchatTEPA, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf("TEPA_final_figures/S09_TEPAsystems_heat.pdf", h = 6, w = 5)
for (i in 1:length(pathways.show.all)) {
netAnalysis_signalingRole_network(cellchatTEPA, signaling = pathways.show.all[i], 
                                  width = 8, height = 2.5, font.size = 10)
}
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("TEPA_project/TEPA-figures/S09_TEPAsystems_dominantSender.pdf", h = 6, w = 5)
gg1 <- netAnalysis_signalingRole_scatter(cellchatTEPA)
gg1
dev.off()
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchatTEPA, signaling = c("CXCL", "CCL"))
gg1 + gg2

pdf("TEPA_project/TEPA-figures/S09_TEPAsystems_dominantAllHeat.pdf", h = 6, w = 12)
ht1 <- netAnalysis_signalingRole_heatmap(cellchatTEPA, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchatTEPA, pattern = "incoming")
ht1 + ht2
dev.off()

# Relate cell groups with their enriched signaling pathways after setting a cutoff for the # pathways for each cell type
# > Cutoff determined by a contribution score
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing") # drops at 9

nPatterns = 8
cellchatTEPA <- identifyCommunicationPatterns(cellchatTEPA, pattern = "outgoing", k = nPatterns)

# river plot
pdf("TEPA_final_figures/S09_TEPAsystems_alluvialOUT.pdf", h = 6, w = 12)
netAnalysis_river(cellchatTEPA, pattern = "outgoing")
dev.off()

selectK(cellchatTEPA, pattern = "incoming")
nPatterns = 8
cellchatTEPA <- identifyCommunicationPatterns(cellchatTEPA, pattern = "incoming", k = nPatterns)
pdf("TEPA_final_figures/S09_TEPAsystems_alluvialIN.pdf", h = 6, w = 12)
netAnalysis_river(cellchatTEPA, pattern = "incoming")
dev.off()

# A - Identify signaling groups based on their functional similarity
cellchatTEPA <- computeNetSimilarity(cellchatTEPA, type = "functional")
cellchatTEPA <- netEmbedding(cellchatTEPA, type = "functional")
cellchatTEPA <- netClustering(cellchatTEPA, type = "functional")
# Visualization in 2D-space
pdf("TEPA_final_figures/S09_TEPAsystems_funct.pdf", h = 6, w = 6)
netVisual_embedding(cellchatTEPA, type = "functional", label.size = 3.5)
dev.off()

# B - Identify signaling groups based on structure similarity
cellchatTEPA <- computeNetSimilarity(cellchatTEPA, type = "structural")
cellchatTEPA <- netEmbedding(cellchatTEPA, type = "structural")
cellchatTEPA <- netClustering(cellchatTEPA, type = "structural")
# Visualization in 2D-space
pdf("TEPA_final_figures/S09_TEPAsystems_struc.pdf", h = 6, w = 6)
netVisual_embedding(cellchatTEPA, type = "structural", label.size = 3.5)
dev.off()

saveRDS(cellchatTEPA, file = "TEPA_results/TEPAcellchat.rds")

