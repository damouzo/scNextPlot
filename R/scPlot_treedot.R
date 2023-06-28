

# Plot TreeDot for GSEA comparision results -----------------------------------
scPlot_treedot <- function(comp_pair_term, top_paths=5, clust_num=3,ORA_type="GO",
                           ORA_ont="BP", ORA_minGSSize=10, ORA_maxGGSSize=500,
                           ORA_p.adj=1, ORA_GO_OrgDb= "org.Hs.eg.db",
                           ORA_p.adj_Meth="BH", ORA_KEGG_Org='hsa'){
  # libraries
  library(ggplot2)
  library(tidyverse)
  library(ggdendro)
  library(cowplot)
  library(clusterProfiler)
  library(ggtree)
  library(patchwork)
  library(org.Hs.eg.db)
  library(RColorBrewer)

  # Estruture Data
  comp_pair_term_fort <- fortify(comp_pair_term, showCategory= top_paths,
                                 includeAll = TRUE, split = NULL)
  comp_pair_term_fort$Cluster <- sub("\n.*", "", comp_pair_term_fort$Cluster)
  comp_pair_term_fort$geneID <- comp_pair_term_fort$core_enrichment

  # Merge Clusters
  merge_compareClusterResult <- function(yy) {
    yy_union <- yy[!duplicated(yy$ID),]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
      ids <- unique(unlist(strsplit(x$geneID, "/")))
      cnt <- length(ids)
      list(ID=paste0(ids, collapse="/"), cnt=cnt)
    })

    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))

    yy_union$geneID <- ids[yy_union$ID]
    yy_union$Count <- cnt[yy_union$ID]
    yy_union$Cluster <- NULL
    yy_union
  }
  merged_ggData <- merge_compareClusterResult(comp_pair_term_fort)

  # Prepare data for IDs Cluster
  prepare_pie_category <- function(enrichDf, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") pie <- "Count"

    pie_data <- enrichDf[,c("Cluster", "Description", "Count")]
    pie_data[,"Description"] <- as.character(pie_data[,"Description"])
    prepare_pie_data(pie_data, pie = pie)
  }
  prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
    if(type == "category"){
      ID_unique <- unique(pie_data[,2])
    } else {
      ID_unique <- unique(pie_data[,3])
    }

    Cluster_unique <- unique(pie_data[,1])
    ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
    if(pie == "Count") {
      for(i in seq_len(nrow(pie_data))) {
        ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
      }
      for(kk in seq_len(ncol(ID_Cluster_mat))) {
        ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
      }
      return(ID_Cluster_mat)
    }
    for(i in seq_len(nrow(pie_data))) {
      if(type == "category"){
        ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
      } else {
        ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
      }

    }
    return(ID_Cluster_mat)
  }
  ID_Cluster_mat <- prepare_pie_category(comp_pair_term_fort,pie = "equal")

  # Hierarchical Clustering
  fill_termsim <- function(x, keep) {
    termsim <- x@termsim[keep, keep]
    termsim[which(is.na(termsim))] <- 0
    termsim2 <- termsim + t(termsim)
    for ( i in seq_len(nrow(termsim2)))
      termsim2[i, i] <- 1
    return(termsim2)
  }
  termsim2 <- fill_termsim(comp_pair_term, rownames(ID_Cluster_mat))
  hc <- stats::hclust(stats::as.dist(1- termsim2),method = "ward.D")

  # Info of Clusters wanted
  clus <- stats::cutree(hc, clust_num)

  # ORA for clusters
  keywords <- vector()
  for (i in 1:clust_num) {
    paths_of_clust <- clus[clus == i]
    cluster_paths_info <- comp_pair_term_fort[comp_pair_term_fort$Description %in% names(paths_of_clust) ,]
    genes2check <-  unique(unlist(str_split(cluster_paths_info$geneID, "/")))

    if (is.null(comp_pair_term@.call$keyType) == F ){
      keytype_called <- comp_pair_term@.call$keyType
    } else if (is.null(comp_pair_term@.call$keyType) == T) {
      paste0("Warning: object@.call$keyType is empty, asumming ENTREZID,
            is and error ocurr please fill keytype argument in compareCluster() ")
      keytype_called <- "ENTREZID"
    }

    # Module for compute ORA base on desired enrichment database
    if (ORA_type == "GO") {
      ora <- enrichGO(genes2check, OrgDb= ORA_GO_OrgDb, keyType=keytype_called,
                      ont= ORA_ont, pvalueCutoff= ORA_p.adj, pAdjustMethod= ORA_p.adj_Meth,
                      qvalueCutoff=1, minGSSize=ORA_minGSSize,maxGSSize=ORA_maxGGSSize)

    } else if (ORA_type == "KEGG") {
      ora <- enrichKEGG(gene= genes2check, organism= ORA_KEGG_Org, pvalueCutoff=ORA_p.adj,
                        pAdjustMethod=ORA_p.adj_Meth, qvalueCutoff=1, minGSSize=ORA_minGSSize,
                        maxGSSize = ORA_maxGGSSize)

    } else if (ORA_type == "DO") {
      ora <- enrichDO(gene=genes2check, ont="DO", pvalueCutoff=ORA_p.adj,
                      pAdjustMethod=ORA_p.adj_Meth, qvalueCutoff=1, minGSSize=ORA_minGSSize,
                      maxGSSize= ORA_maxGGSSize)
    }

    keywords[i] <- ora@result$Description[1]
  }

  #For leyend adjust proportionaly to terms length
  keywords_max_length =max(nchar(as.list(keywords)))
  if (keywords_max_length < 17){keywords_max_length = 17}

  # Plot dendogram
  g <- split(names(clus), clus)
  p <- ggtree(hc, size=1.2)
  clades <- sapply(g, function(n) MRCA(p, n))
  p <- groupClade(p, clades, group_name='SubTree_ORA') + aes(color=SubTree_ORA)
  mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 6))
  ggtree_plot <- p + scale_color_manual(values=mycolors, breaks=1:clust_num,labels = keywords)+
    theme(legend.position='right',legend.justification = c(0,1.5))
  ggtree_plot_noLegend <- ggtree_plot + theme(legend.position = "none")


  # Prepare for dotplot
  comp_pair_term_fort$log_p.adjust <- -log10(comp_pair_term_fort$p.adjust)
  comp_pair_term_fort$Description = factor(comp_pair_term_fort$Description,
                                           levels = hc$labels[hc$order])
  comp_pair_term_fort$Cluster <- factor(comp_pair_term_fort$Cluster,
                                        levels=levels(comp_pair_term@compareClusterResult$Cluster))

  # Plot dotplot
  dotplot <- comp_pair_term_fort %>%
    ggplot(aes(x=Cluster, y =Description , color = NES, size = log_p.adjust)) +
    geom_point() +
    scale_y_discrete(position = "right")+
    scale_color_gradient2(low="blue4",mid="white", high="red")+
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') +
    guides(size=guide_legend(title="-log10(p.adj)"))+
    theme(axis.ticks = element_blank(),legend.position = "right",legend.justification = c(0,0))
  dotplot_noLegend <- dotplot + theme(legend.position = "none")

  # Legends
  legend_tree <- get_legend(ggtree_plot)
  legend_dot <- get_legend(dotplot)

  # Merge Plot
  combine_plot <- plot_grid(ggtree_plot_noLegend,NULL, dotplot_noLegend, nrow= 1, rel_widths= c(0.3,-0.05,2), align = 'h')
  combine_legend <- plot_grid(legend_dot,NULL,legend_tree, ncol=1, rel_heights = c(1,-0.5,1))
  big_plot <- plot_grid(combine_plot,combine_legend,NULL, nrow = 1,
                        rel_widths = c(1,0.1, keywords_max_length*0.006))
  return(big_plot)
}

