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
library("CellChat") # Run with v2.0
library(patchwork)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 

setwd("~/Library/CloudStorage/OneDrive-UNSW/TEPA_project")
source("TEPA_code/supportFunctions.R")
#seuset_con <- LoadSeuratRds("TEPA_results/S08_seusetFull.Rds")

# 1. DATA PREPARATION ####

seuset_con <- seuset_con[,seuset_con$condition == "Control"]

data.input <- GetAssayData(seuset_con, assay = "RNA", layer = "data") # normalized data matrix
labels <- seuset_con$celltypes
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchatCON <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchatCON <- addMeta(cellchatCON, meta = meta)
cellchatCON <- setIdent(cellchatCON, ident.use = "group") # set "labels" as default cell identity
levels(cellchatCON@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchatCON@idents)) # number of cells in each cell group

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
cellchatCON@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchatCON <- subsetData(cellchatCON) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchatCON <- identifyOverExpressedGenes(cellchatCON)
cellchatCON <- identifyOverExpressedInteractions(cellchatCON)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchatCON <- projectData(cellchatCON, PPI.human)

# 3. INFERENCE OF CELL-CELL COMMUNICATION METHOD ####

# Compute the communication probability and infer cellular communication network
cellchatCON <- computeCommunProb(cellchatCON, type = "triMean", population.size = TRUE)
#> triMean is used for calculating the average gene expression per cell group. 
#> ‘trimean’ approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. >> most important interactions very stringent
#> To use 10% truncated mean, USER can set type = "truncatedMean" and trim = 0.1. 
#> To determine a proper value of trim, function computeAveExpr can check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchatCON, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1)

cellchatCON <- filterCommunication(cellchatCON, min.cells = 10)

df.net <- subsetCommunication(cellchatCON)
df.net_short <- subsetCommunication(cellchatCON, sources.use = c("Neutrophils"), targets.use = c("Gamma-delta T cells","Cd4+ Naive T cells", "DN Regulatory T cells",
                                                                                              "Cd8+ Naive-Memory T cells", "Cd8+ NkT-like cells", "Natural killer cells"))

#> df.net <- subsetCommunication(cellchatCON, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

cellchatCON <- computeCommunProbPathway(cellchatCON)

cellchatCON <- aggregateNet(cellchatCON)

groupSize <- as.numeric(table(cellchatCON@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("TEPA_project/TEPA-figures/S09_CONcellComm1.pdf", h = 10, w = 10)
netVisual_circle(cellchatCON@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatCON@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("TEPA_project/TEPA-figures/S09_CONcellComm2.pdf", h = 10, w = 15)
mat <- cellchatCON@net$weight
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
pdf("TEPA_final_figures/S09_CONcellComm_TNF.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_aggregate(cellchatCON, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# B- Chord diagram -> does it take into account population proportions?
par(mfrow=c(1,1))
pathways.show <- c("CXCL") 
pdf("TEPA_final_figures/S09_CONcellComm_CXCL_chord.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_aggregate(cellchatCON, signaling = pathways.show, layout = "chord")
dev.off()

par(mfrow=c(1,1))
pathways.show <- c("TNF") 
pdf("TEPA_final_figures/S09_CONcellComm_TNF_chord.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_aggregate(cellchatCON, signaling = pathways.show, layout = "chord")
dev.off()

# C - Heatmap
par(mfrow=c(1,1))
pathways.show <- c("CXCL") 
pdf("TEPA_final_figures/S09_CONcellComm_TNF_heatmap.pdf", h = 10, w = 10) # TNF signalling pathway network
netVisual_heatmap(cellchatCON, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

# see contribution of each LR pair
pathways.show <- c("TNF")  # Tnfrsf1b/a
netAnalysis_contribution(cellchatCON, signaling = pathways.show) # you can include only one LR pair in the analysis

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchatCON@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchatCON@idents)
vertex.receiver = seq(1,4)
pdf("TEPA_final_figures/S09_CONcellComm_All_chord.pdf", h = 10, w = 10)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_aggregate(cellchatCON, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "chord")
}
dev.off()

pdf("TEPA_final_figures/S09_CONcellComm_All_circle.pdf", h = 12, w = 10)
for (i in 1:length(pathways.show.all)) {
  netVisual_aggregate(cellchatCON, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver,
                      signaling.name = pathways.show.all[i])
}
dev.off()

# Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
for (i in 1:length(pathways.show.all)) {
  p<-netAnalysis_contribution(cellchatCON, signaling = pathways.show.all[i])
  file = paste0("TEPA_final_figures/S09_CONcellComm_contribution_", pathways.show.all[i])
  ggsave(p, filename = file, width = 4, height = 2, units = 'in', dpi = 300, device = "pdf")
}

# D- Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
levels(cellchatCON@meta$group)
pdf("TEPA_project/TEPA-figures/S09_CONcellComm_neutrophils_vs_all_BUBBLE.pdf", h = 6, w = 5)
netVisual_bubble(cellchatCON, sources.use = "Neutrophils", targets.use = c(1:14), remove.isolate = FALSE)
dev.off()

pdf("TEPA_project/TEPA-figures/S09_CONcellComm_all_vs_neutrophils_BUBBLE.pdf", h = 6, w = 5)
netVisual_bubble(cellchatCON, sources.use = c(1:14), targets.use = "Neutrophils", remove.isolate = FALSE)
dev.off()

# 5. Systems analysis of cell-cell communication network ####

# Compute the network centrality scores
cellchatCON <- netAnalysis_computeCentrality(cellchatCON, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf("TEPA_final_figures/S09_CONsystems_heat.pdf", h = 6, w = 5)
for (i in 1:length(pathways.show.all)) {
netAnalysis_signalingRole_network(cellchatCON, signaling = pathways.show.all[i], 
                                  width = 8, height = 2.5, font.size = 10)
}
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("TEPA_project/TEPA-figures/S09_CONsystems_dominantSender.pdf", h = 6, w = 5)
gg1 <- netAnalysis_signalingRole_scatter(cellchatCON)
gg1
dev.off()
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchatCON, signaling = c("CXCL", "CCL"))
gg1 + gg2

pdf("TEPA_project/TEPA-figures/S09_CONsystems_dominantAllHeat.pdf", h = 6, w = 12)
ht1 <- netAnalysis_signalingRole_heatmap(cellchatCON, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchatCON, pattern = "incoming")
ht1 + ht2
dev.off()

# Relate cell groups with their enriched signaling pathways after setting a cutoff for the # pathways for each cell type
# > Cutoff determined by a contribution score
library(NMF)
library(ggalluvial)
selectK(cellchatCON, pattern = "outgoing") # drops at 9

nPatterns = 8
cellchatCON <- identifyCommunicationPatterns(cellchatCON, pattern = "outgoing", k = nPatterns)

# river plot
pdf("TEPA_final_figures/S09_CONsystems_alluvialOUT.pdf", h = 6, w = 12)
netAnalysis_river(cellchatCON, pattern = "outgoing")
dev.off()

selectK(cellchatCON, pattern = "incoming")
nPatterns = 8
cellchatCON <- identifyCommunicationPatterns(cellchatCON, pattern = "incoming", k = nPatterns)
pdf("TEPA_final_figures/S09_CONsystems_alluvialIN.pdf", h = 6, w = 12)
netAnalysis_river(cellchatCON, pattern = "incoming")
dev.off()

# A - Identify signaling groups based on their functional similarity
cellchatCON <- computeNetSimilarity(cellchatCON, type = "functional")
cellchatCON <- netEmbedding(cellchatCON, type = "functional")
cellchatCON <- netClustering(cellchatCON, type = "functional")
# Visualization in 2D-space
pdf("TEPA_final_figures/S09_CONsystems_funct.pdf", h = 6, w = 6)
netVisual_embedding(cellchatCON, type = "functional", label.size = 3.5)
dev.off()

# B - Identify signaling groups based on structure similarity
cellchatCON <- computeNetSimilarity(cellchatCON, type = "structural")
cellchatCON <- netEmbedding(cellchatCON, type = "structural")
cellchatCON <- netClustering(cellchatCON, type = "structural")
# Visualization in 2D-space
pdf("TEPA_final_figures/S09_CONsystems_struc.pdf", h = 6, w = 6)
netVisual_embedding(cellchatCON, type = "structural", label.size = 3.5)
dev.off()

saveRDS(cellchatCON, file = "TEPA_results/TEPAcellchatCON.rds")

