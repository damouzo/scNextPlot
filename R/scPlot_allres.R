
# Compile descriptive resolutions ----------------------------------------------
scPlot_allres <- function(Seu_obj,res=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                          clustree=T, UMAPs=T, clustreeGenes=NULL ,wknnGraph=F) {
  # Load required libraries
  library(Seurat)
  library(clustree)
  library(ggplot2)
  library(dplyr)
  library(BuenColors)

  #Load Colors
  my_pal <- c(jdb_palette("corona"),jdb_palette("lawhoops"),
              jdb_palette("samba_color"),jdb_palette("wolfgang_basic"))

  # Clean Res and Compute FindCluster for selected resolutions ------------------
  Seu_obj@meta.data <- Seu_obj@meta.data %>% dplyr::select(-contains("res"))

  if (wknnGraph == T) {
    Seu_obj <- FindClusters(Seu_obj, resolution = res, graph="wknn")

  } else {
    Seu_obj <- FindClusters(Seu_obj, resolution = res)
  }


  #Create list of plots for disfferent resolutions -----------------------------
  col_names <- grep("res", names(Seu_obj@meta.data), value = TRUE)
  plots <- vector('list', length=length(col_names))
  names(plots) <- col_names
  # Create a DimPlot for each resolution
  if (UMAPs == T){
    for (plot in names(plots)) {
      plots[[plot]] <- DimPlot(Seu_obj, reduction = "umap", group.by = plot,
                               label = T,cols= my_pal)
    }
  }

  # Create a clustree for summary resolution flow ------------------------------
  # Set config
  prefix <- unique(sub("\\..*", ".", col_names))
  P1 <- NULL
  plots_cgenes <- NULL

  if (clustree == TRUE){
    P1 <- clustree(Seu_obj@meta.data, prefix = prefix) # Print Plot
  }

  # Create clustree an check specific genes ------------------------------------
  if (!is.null(clustreeGenes)){
    # make List for plot for each gene
    plots_cgenes <- vector('list', length=length(clustreeGenes))
    names(plots_cgenes) <- clustreeGenes

    #Get Gene Info and add to metadata
    genes_res<- t(as.data.frame(Seu_obj@assays$RNA@data[rownames(Seu_obj@assays$RNA@data) %in% clustreeGenes , , drop=F]))
    Seu_obj@meta.data <- merge(Seu_obj@meta.data,genes_res,by='row.names')

    for (cgene in colnames(genes_res)) {
      plots_cgenes[[cgene]] <- clustree(Seu_obj@meta.data, prefix= prefix,
                                        node_colour= cgene,node_colour_aggr = "mean") +
        scale_colour_viridis_c(option = 'plasma', begin = 0.3) +
        ggtitle(cgene)
    }
  }

  # Plot if exist
  if (!is.null(P1)) {print(P1)}
  if (!is.null(plots_cgenes[[1]])) {print(plots_cgenes)}
  if (!is.null(plots[[1]])) {print(plots)}

}

