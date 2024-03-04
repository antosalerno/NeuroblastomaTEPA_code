# Define a set of supporting function that will be used in the pipeline #
## author: Antonietta Salerno
## date: 20/01/2023

# Define color coding scheme
cellt_col <- c("#CC4C02","#CE1256","#004529","#ABA9E8",
                        "#5AAE61","#C979D7", "#6BAED6" , 
                         "#FED976","#7A0177","#35978F","#67001F",
                        "#045A8D", "#A6C950","#B30000")
                        
tum_col <- colorRampPalette(RColorBrewer::brewer.pal(10,"Paired"))(10)[3:6]
cond_col <- c("#ADD8E6", "#A50F15")
cond_col <- colorRampPalette(RColorBrewer::brewer.pal(4,"Paired"))(4)[1:2]
names(cond_col) <- c("Control", "Treatment")

inf_col <- colorRampPalette(RColorBrewer::brewer.pal(4,"Paired"))(4)[3:4]
names(inf_col) <- c("T", "F")

# S00 - Density histogram ####

hist_dens <- function(x, breaks = "Scott", 
                      main = "Mycn expression per number of cells",
                      xlab = "Mycn", ylab = "Cell counts") {
  dens <- density(x, na.rm = T)
  raw_hist <- hist(x, breaks = breaks, plot = F)
  scale <- max(raw_hist$counts)/max(raw_hist$density)
  hist(x, breaks = breaks, prob = F, main = main, xlab = xlab, ylab = ylab)
  lines(list(x = dens$x, y = scale * dens$y), col = "blue", lwd = 2)
}

# S03 - Differential expression analysis ####

# Create signature of cell type with top markers
createSets <- function(markers = immune.markers,
                       obj = seuset_immune,
                       id = "celltypes"){
  Idents(obj) <- id
  clusters = unique(Idents(obj))
  for(c in 1:length(clusters)){
    clusterDF <- markers[markers$cluster == clusters[c],]
    if (nrow(clusterDF) > 0){
      clusterDF <- clusterDF[clusterDF$p_val_adj < 0.05,]
      clusterDF <- clusterDF[order(-c(clusterDF$avg_log2FC)),]
      top_genes <- head(clusterDF[clusterDF$cluster == clusters[c],]$gene)
      obj <- AddModuleScore(obj, assay = "RNA", features = list(top_genes), name=make.names(as.character(clusters[c])))
      names(obj@meta.data)[grep(make.names(as.character(clusters[c])), names(obj@meta.data))] <- as.character(clusters[c])
    }
  }
  return(obj)
}

# Get top marker genes by cluster
getTopMarkers <- function(df = immune.markers, geneNr){
  markers <- c()
  for(c in 1:length(clusters)){
    clusterDF = df[df$cluster == clusters[c],]
    clusterDF <- clusterDF[order(-clusterDF$avg_log2FC),]
    top_genes <- head(clusterDF[clusterDF$cluster == clusters[c],]$gene, geneNr)
    markers <- append(markers, top_genes)
  }
  return(markers)
}

# Volcano plot

library("EnhancedVolcano")

plotVolcano <- function(clusters, res=NULL, type = "condition",
                        save = "", log2FC = 0.25, immune = TRUE, pval = 0.05){
  for (cluster in clusters){
    try({
      if (type == "condition"){
        if(immune){
          file=paste0("TEPA_results/S03_immuneCond_DEA_",gsub(" |/", "_", cluster),".csv")
        } else {
          file=paste0("TEPA_results/S05_tumorCond_DEA_",gsub(" |/", "_", cluster),".csv")
        } 
        res <- read.csv(file, sep=",")
        rownames(res) <- res$X
        title = paste0(cluster,', TEPA vs. CTRL ')
      } else if (type == "markers") {
        if (immune){res <- immune.markers[immune.markers$cluster == cluster,]}
        else {res <- tumor.markers[tumor.markers$cluster == cluster,]}
        title = paste0('DEA cluster ', cluster, " vs all")
      } else {message("Please specify the type of analysis: either 'markers' or 'condition'")}
      p <- EnhancedVolcano(res, subtitle = "",
                         #selectLab = markers,
                         lab = rownames(res),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         xlim = c(-2.5, 2.5),
                         title = title,
                         pCutoff = pval, 
                         FCcutoff = log2FC, # 2-fold change
                         labFace = "bold",
                         labSize = 3,
                         col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                         colAlpha = 4/5,
                         legendLabSize = 14,
                         legendIconSize = 4.0,
                         drawConnectors = TRUE,
                         widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                         caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC > log2FC&res$p_val_adj<= 0.05,]), ' genes',
                                          "\n",'Downregulated = ', nrow(res[res$avg_log2FC < -log2FC &res$p_val_adj<= 0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
      ggsave(p, file=paste0("TEPA_plots/",save, "DEA_",sub(" ", "_", cluster),".png"), width = 30, height = 25, units = "cm")
    })
  }
}

# S04 - GSEA ####

N1 <- c("Cebpb", "Il1b","Cxcl2","S100a8","S100a9","Ccl3", "Ccl4", "Ifitm1", "Ifitm2", "Ifitm3", "Ifitm6", 
          "Acod1","Sell", "Prkcd","Il6ra","Retnlg","Il1r2","Mmp8", "Mmp9","Icam1","Hif1a",  
         "Tnf", "Myd88", "Fas", "Cxcl3", "Isg15", "Isg20", "Arg2", "Stat3", 
         "Il1a" , "Il15")

N2 <- c("Tgfb1", "Tgfb1i1","Tgfb2", "Tgfb3", "Ccl2",  "Ccl17","Cxcl14", 
        "Cxcl15",  "Il1r1", "Il2", "Il17a", "Mpo", "Slc27a2", "Arg1", 
        "Mrc1", "Chil3", "Elane", "Ctsg", "Retnla", "Siglecf")

# https://www.sciencedirect.com/science/article/pii/S1535610809002153?via%3Dihub
# https://www.sciencedirect.com/science/article/pii/S1471490605002425
# https://ashpublications.org/blood/article/133/20/2159/273823/Neutrophil-plasticity-in-the-tumor
# https://urldefense.com/v3/__https://pubmed.ncbi.nlm.nih.gov/37001504/__;!!I9nncg!r5yciw1EwGP_GPPv2xDM4awc2iWXGyZrOLYVfYyun97DXla3h8RIigjJfJr93_YbCYZZUG74WDyzVrFnNwnB$

M1 <- c("Il12a", "Il12b", "Tnf", "Cxcl10", "Tlr2", "Tlr4", "Fcgr3", "Fcgr4", "Fcgr1",
        "Cd80", "Cd86", "Il12a",  "Il6", "Il1a", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
        "Cxcl9", "Cxcl10", "Ccr7", "Il1b", "Nos2", "Cxcl16", "Il23a")
M2 <- c("Il10", "Ccl22", "Il4ra", "Il13ra1", "Chil3", "Cd163", "Tgfb1",
        "Fcer2a", "Cd163", "Cxcr1", "Cxcr2", "Ccr2", "Arg1",  "Cd209a", "Mrc1", "Fcgr2b")

# https://www.sciencedirect.com/science/article/pii/S1471490602023025

custom <- list(N1, N2, M1, M2)
names(custom) <- c("N1 ANTI-TUMOR NEUTROPHILS", "N2 PRO-TUMOR NEUTROPHILS",
                   "M1 ANTI-TUMOR MACROPHAGES", "M2 PRO-TUMOR MACROPHAGES")

copper_genes <- c("Slc31a1", "Atp7a", "Atp7b", "Sco1", "Cox11", "Steap3", "Commd1", "Mtf1", "Mtf2", "Sp1", "Sod1",
                  "Sod2", "Steap4", "Atox1", "Ccs", "Mt1", "Mt2", "Mt3", "Fdx1", "Lias", "Lipt1", "Dld", "Dlat",
                  "Pdha1", "Pdhb",  "Gls", "Cdkn2a")

copper_sign <- list(copper_genes)
names(copper_sign) <- "COPPER METABOLISM"

# @ version of clusterProfiler:: cnetplot
networkPlotGSEA <- function(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 6){
  out <- tryCatch({
    cnet <- clusterProfiler::cnetplot(fgseaResPlot, showCategory = CNET_gene_cutoff, colorEdge = TRUE, 
                                      categorySize="p.adjust", color.params = list(foldChange = ranked.genes)) + 
      labs(title = paste("Significant Pathways in ", cluster ,sep = ""),
           subtitle = "enriched by GSEA",
           caption = paste0("Based on ",
                            length(ranked.genes[ranked.genes > 0])," upregulated genes and ",
                            length(ranked.genes[ranked.genes < 0])," downregulated genes.",
                            " Showing ",CNET_gene_cutoff," genes per pathway."))+
      theme(plot.title = element_text(hjust = 0.5, size=20),
            plot.subtitle = element_text(hjust = 0.5, size=15),
            plot.caption = element_text(size=12))
    },
    finally = {
      return(cnet)
    })
  return(out)
}

library(scales)

barPlotGSEA <- function(fgseaRes, cluster = NULL, name = "", byType = FALSE){
  
  # When you want to plot the enrichment of gene sets by cell type
  if (byType){
    
    #clusters_colors = hue_pal()(length(clusters))
    #order_cl <- sort(unlist(clusters))
    #c <- setNames(clusters_colors, order_cl)
    c = cellt_col[1:13]
    
    fgseaRes %>% arrange(pathway, NES) %>%
      ggplot(aes(x = pathway, 
                 y = as.numeric(NES),
                 fill = reorder(cluster, as.numeric(NES)))) +
      scale_colour_manual(values = c) +
      geom_bar(position=position_dodge2(width = 1.5, preserve = "single", padding = 0),
               stat='identity', width = 1) +
      ylab("Normalized enrichment score") +
      geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
      scale_x_discrete(expand = c(.05, .05)) +
      coord_cartesian(xlim = c(-5, 5)) +
      coord_flip() +
      # scale_y_discrete(limits = factor(min(fgseaRes$NES), max(fgseaRes$NES)), breaks = seq(min(fgseaRes$NES), max(fgseaRes$NES),by = 1)) +
      labs(title = "Top pathways enriched by cell type",
           x = NULL,
           fill = NULL) +
      theme_minimal() +
      theme(panel.grid.major.y = element_blank(),
            legend.position = "right") +
      theme(axis.text.y = element_text(size = 10, hjust = 1),
            plot.margin = margin(rep(10, 8)))
  }
  
  # When you want to plot the enrichment of a specific cluster
  #if (is.null(cluster) == TRUE)
  else {
  diverging_bar_chart(fgseaRes, pathway, as.numeric(NES), 
                      bar_color = c("#F9AF93", "#9BD6F0"),
                      text_color = "black", text_size = 10) +
    ylab("Normalized enrichment score") +
    geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
    ggtitle(paste0("Significant pathways in ", name," enriched by GSEA")) + 
    theme(plot.title = element_text(hjust = 1)) + 
    theme(axis.title=element_text(size=20)) +
    theme_minimal() +
    ylim(-5, 5) +
    theme(axis.text.y = element_blank(),  # Remove Y-axis texts
          axis.ticks.y = element_blank(), # Remove Y-axis ticks
          panel.grid.major.y = element_blank()) # Remove horizontal grid
  }
}


# Run fgsea function from DEA results using GO database and produce networkPlot and barPlot 
gseaRES <- function(clusters, markers = NULL, fgsea_sets, minSize = 12, adj = TRUE, 
                    type = "condition", input = "immune", save="", out = "png"){
  for (cluster in clusters){
    try({
      if (type == "condition"){
        if(input == "immune"){
          #file=paste0("TEPA_results/S03_immuneBulkDEA.csv")
          file=paste0("TEPA_results/S03_immuneCond_DEA_",gsub(" ", "_", cluster),".csv")
        } else if (input == "tumor"){
          file=paste0("TEPA_results/S05_tumorBulkDEA.csv")
          #file=paste0("TEPA_results/S05_tumorCond_DEA_",gsub(" |/", "_", cluster),".csv")
        } else if (input == "nanoCond_gInf"){
            file = paste0("TEPA_results/N00_nanoCond_gInf_DEA.csv")
        } else if (input == "nanoInf_gCond"){
            file = paste0("TEPA_results/N00_nanoInf_gCond_DEA.csv")
        } else if (input == "nanoCond"){
           file = paste0("TEPA_results/N00_nanoCond_DEA.csv")
        }
        res <- read.csv(file = file, sep=",")
      } else if (type == "markers") {
        res <- markers[markers$cluster == cluster,]
        res$X <- res$gene
      } else {message("Please specify the type of analysis: either 'markers' or 'condition'")}

      message(paste0("GSEA for cluster: ", cluster))
      ranked.genes.df<- res %>%
                        arrange(desc(avg_log2FC)) %>% 
                        dplyr::select(avg_log2FC, X)
      ranked.genes <- ranked.genes.df$avg_log2FC
      names(ranked.genes) <- ranked.genes.df$X
      ranked.genes <- ranked.genes[!is.na(names(ranked.genes))]
    
      # Run GSEA
      fgseaRes<- fgsea(fgsea_sets, stats = ranked.genes, minSize = minSize, maxSize = Inf)
      if (adj) {fgseaRes$padj = p.adjust(fgseaRes$pval, method='BH')}
      fgseaRes <- fgseaRes[fgseaRes$padj <= 0.1] %>% arrange(desc(NES))
      fgseaRes <- as.data.frame(fgseaRes)
      
      # enr <- plotEnrichment(fgsea_sets[["N1 ANTI-TUMOR NEUTROPHILS"]],
      #                ranked.genes) + labs(title="N1 ANTI-TUMOR NEUTROPHILS") + geom_line(color="#045A8D")
      # #NES = -1.7
      # 
      # file=paste0("TEPA_final_figures/S05_N1_NEUTROPHILS.pdf")
      # ggsave(enr, file = file, width = 10, height = 10, units = "cm")
      
    
      # Generate enrichment plots 
      if (nrow(fgseaRes) > 0) {
        fgseaResPlot <- fgseaRes %>% split(x = .$leadingEdge, 
                                         f = .$pathway)
        for (e in 1: length(fgseaResPlot)){
          fgseaResPlot[[e]] <- unlist(fgseaResPlot[[e]])
        }
      
      ## Network
      cnet <- suppressWarnings(networkPlotGSEA(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 10))
      file=paste0("TEPA_plots/",save, "cnet_",sub(" ", "_", cluster),".", out)
      ggsave(cnet, file = file, width = 30, height = 30, units = "cm")
      
      ## Barplot
      b <- barPlotGSEA(fgseaRes, cluster, name = cluster)
      file=paste0("TEPA_plots/",save, "barplot_",sub(" ", "_", cluster),".", out)
      ggsave(b, file = file,width = 20, height = ifelse(nrow(fgseaRes) > 4, nrow(fgseaRes)/1.5, nrow(fgseaRes) + 4), units = "cm")
      
      # Save dataframe
      df <- apply(fgseaRes,2,as.character)
      file=paste0("TEPA_results/",save, "GSEA_",sub(" ", "_", cluster),".csv")
      write.csv(df, file = file)
      }
    }
  )}
}

gseaPlotRes <- function(clusters, save="", out="png"){
  for (cluster in clusters){
    message(paste0("GSEA for cluster: ", cluster))
    
    # Load DEA analysis
    res <- read.csv(paste0("TEPA_results/",save, cluster,".csv"), sep=",")
    ranked.genes.df<- res %>%
      arrange(desc(avg_log2FC)) %>% 
      dplyr::select(avg_log2FC, X)
    ranked.genes <- ranked.genes.df$avg_log2FC
    names(ranked.genes) <- ranked.genes.df$X
    ranked.genes <- ranked.genes[!is.na(names(ranked.genes))]
    
    # Load GSEA analysis
    filename = paste0("TEPA_results/",save, cluster,"SHORT.csv")
    
    if (file.exists(filename)){
      fgseaRes <- read.csv(filename, sep=";")
    
      fgseaResPlot <- fgseaRes %>% split(x = .$leadingEdge, 
                                         f = .$pathway)
      for (e in 1: length(fgseaResPlot)){
        leEd <- unlist(str_extract_all(fgseaRes$leadingEdge[e],"[[:alnum:]]{2,}"))
        fgseaResPlot[e][[1]] <- leEd
      }
      
      ## Network
      cnet <- suppressWarnings(networkPlotGSEA(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 10))
      ggsave(cnet, file=paste0("TEPA_plots/",paste(cluster, collapse = "_"),save,"_cnet.", out), width = 50, height = 50, units = "cm")
      
      ## Barplot
      b <- barPlotGSEA(fgseaRes, cluster)
      ggsave(b, file=paste0("TEPA_plots/",paste(cluster, collapse = "_"),save,"_barplot.", out),width = 50, height = ifelse(nrow(fgseaRes) > 3, nrow(fgseaRes)/1.5, nrow(fgseaRes) + 1), units = "cm")
    }
  }
}
  

gseaByType <- function(clusters, save, custom_sign = FALSE){
  
  fgseaResByType <- data.frame(
    pathway = character(0),
    NES = numeric(0),
    cluster = character(0)
  )
  for (cluster in clusters){
    index = which(clusters == cluster)
    if (custom_sign){
      filename = paste0("TEPA_results/S04_immuneEnrichmentCustom_GSEA_",sub(" ", "_", cluster),".csv")
    } else{
      filename = paste0("TEPA_results/S04_immuneEnrichment_GSEA_",sub(" ", "_", cluster),".csv")
    } 
    if (file.exists(filename)){
      clusterFile <- read.csv(filename, sep = ",")
      if (ncol(clusterFile) == 2){
        fgseaResByType[nrow(fgseaResByType)+1,"pathway"] <- clusterFile[1,2]
        fgseaResByType[nrow(fgseaResByType),"NES"] <- clusterFile[6,2]
        fgseaResByType[nrow(fgseaResByType),"cluster"] <- cluster
      }
      else {
        for (i in 1: nrow(clusterFile)){
          fgseaResByType[nrow(fgseaResByType)+1,]$pathway <- clusterFile[i,"pathway"]
          fgseaResByType[nrow(fgseaResByType),]$NES <- clusterFile[i,"NES"]
          fgseaResByType[nrow(fgseaResByType),]$cluster <- cluster
        }
      }
    }
  }
  
  # Save dataframe
  #df <- apply(fgseaResByType,2,as.character)
  write.csv(fgseaResByType, file = paste0("TEPA_results/", save, ".csv"))
  
}

 # Lots of time -> try Reactome as well
gseaJointNet <- function (clusters, immune = TRUE, save=""){
  cell_types = list()
  for (cluster in clusters){
    if(immune){
      file=paste0("TEPA_results/S03_immuneCond_DEA_",gsub(" ", "_", cluster),".csv")
    } else {
      file=paste0("TEPA_results/S07_tumorCond_DEA_",gsub(" ", "_", cluster),".csv")
    } 
    res <- read.csv(file = file, sep=",")
    ranked.genes.df<- res %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(avg_log2FC > 0.2) %>%
      arrange(desc(avg_log2FC)) %>% 
      dplyr::select(avg_log2FC, X)
    ranked.genes <- ranked.genes.df$avg_log2FC
    names(ranked.genes) <- ranked.genes.df$X
    
    # We take only significant positively enriched genes since we are running ORA
    if (length(ranked.genes) > 0){
      cell_types[cluster] <- list(names(ranked.genes))
    
    }
  }
  
  gseaByCellType <- compareCluster(cell_types, fun="enrichGO", ont = "BP", readable = TRUE, keyType = "SYMBOL",
                                   OrgDb = "org.Mm.eg.db", pvalueCutoff=0.05, minGSSize = 15, maxGSSize = 500)
  
  df_gseaCT <- gseaByCellType@compareClusterResult
  df_gseaCT$geneNum <- as.numeric(str_extract_all(df_gseaCT$GeneRatio,"[[:digit:]]{1,}", simplify = T)[,1])
  df_gseaCT <- df_gseaCT %>%
    filter(p.adjust <=0.05) %>%
    filter(geneNum > 10) # remove gene sets matches having less than 10 genes 
  
  df_gseaCT <- as.data.frame(df_gseaCT)
  
  # Save dataframe
  df <- apply(df_gseaCT,2,as.character)
  file=paste0("TEPA_results/",save, ".csv")
  write.csv(df, file = file)
  
  return(gseaByCellType)
}

# Plotting gene sets ####
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

sign_dotPlot <- function(obj, sign, immune = TRUE, legend = TRUE, k = 3, cluster = TRUE){
  Idents(obj) <- "celltypes"
  p <- DotPlot(object = obj, features = sign, 
               scale=FALSE, col.min = -4, col.max = 4, assay = "RNA",
               dot.min = 0, dot.scale = 2, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
    theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
  df <- p$data
  
  ### the matrix for the scaled expression 
  exp_mat<-df %>% 
    select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame() 
  
  row.names(exp_mat) <- exp_mat$features.plot  
  exp_mat <- exp_mat[,-1] %>% as.matrix()
  
  percent_mat<-df %>% 
    select(-avg.exp, -avg.exp.scaled) %>%  
    pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame() 
  
  row.names(percent_mat) <- percent_mat$features.plot  
  percent_mat <- percent_mat[,-1] %>% as.matrix()
  
  quantile(exp_mat, c(0.1, 0.3, 0.7, 0.99))
  
  cell_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= (percent_mat[i, j]/100 * min(unit.c(w, h)))/1.5,
                gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
  names(cellt_col) <- levels(obj$celltypes)
  
  if (immune){
    ha = HeatmapAnnotation(show_legend = legend,
                           celltypes = levels(obj$celltypes),
                           col = list("celltypes" = cellt_col[1:13]))
  } else {
    ha = HeatmapAnnotation(show_legend = legend,
                           celltypes = levels(obj$celltypes),
                           col = list("celltypes" = cellt_col))
  }
  
  col_fun = colorRamp2(seq(min(exp_mat), max(exp_mat), length.out = 2), c("#EEEEEE", "red"))
  Heatmap(exp_mat,
          heatmap_legend_param=list(title="z-score"),
          col=col_fun,
          #clustering_distance_rows = "pearson",
          rect_gp = gpar(type = "none"),
          cell_fun = cell_fun,
          cluster_columns = FALSE,
          cluster_rows = cluster,
          column_order = NULL,
          row_names_gp = gpar(fontsize = 15, color = "white", lwd = 2),
          row_title = NULL,
          row_km = ifelse(cluster,k,1),
          use_raster = TRUE,
          raster_quality = 5,
          raster_resize_mat = TRUE,
          column_names_rot = 45,
          row_dend_reorder = TRUE,
          top_annotation = ha,
          heatmap_width = unit(25, "cm"), 
          heatmap_height = unit(30, "cm"),
          border_gp = gpar(col = "black", unit(10, "mm")))
  
}


sign_avgHeatMap <- function(obj, sign, immune = TRUE, 
                            legend = TRUE, k = 3, cluster = TRUE,
                            w = 25, h = 30){
  
  # Expression data
  my_data <- AverageExpression(
    obj,
    assays = "RNA",
    features = sign,
    group.by = c("celltypes", "condition"),
    layer = "data")$RNA
  
  my_data_matrix = as.matrix(my_data)
  rownames(my_data_matrix) = rownames(my_data)
  colnames(my_data_matrix) = colnames(my_data)
  
  # order of annotations/colors are defined here
  ordered_meta_data <- str_split_fixed(colnames(my_data), '_', 2)
  rownames(ordered_meta_data) <- colnames(my_data)   
  colnames(ordered_meta_data) <- c("celltypes", "condition")
  
  #levels(cellt_col) <- levels(seuset_full$celltypes)
  names(cellt_col) <- levels(obj$celltypes)
  
  if (immune){
    annotation_colors <- list("celltypes" = cellt_col[1:13],
                              "condition" = cond_col)
  } else{
    annotation_colors <- list("celltypes" = cellt_col,
                              "condition" = cond_col)
  }
  
  ha = HeatmapAnnotation(df = as.data.frame(ordered_meta_data),
                         show_legend = legend,
                         show_annotation_name = TRUE,
                         col = annotation_colors,
                         annotation_legend_param = list(
                           celltypes = list(
                             labels = levels(obj$celltypes)),
                           title = "celltypes"
                         ))
  
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  Heatmap(
    my_data_matrix,
    col = col_fun,
    cluster_rows = cluster,
    row_km = ifelse(cluster,k,1),
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
    heatmap_width = unit(w, "cm"), 
    heatmap_height = unit(h, "cm")
    #raster_by_magick = TRUE,
    #raster_quality = 5,
    #raster_resize_mat = TRUE
  )
  
}

