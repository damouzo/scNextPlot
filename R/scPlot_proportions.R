
# Proportions ident in scRNA plot ---------------------------------------------
scPlot_proportions = function (seurat_obj, by_ident, clusters=F,
                               PlotNcells=T, value_reorder=T) {
  ## take an Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(BuenColors)

  count_table <- table(seurat_obj@meta.data[[by_ident]], seurat_obj@meta.data$orig.ident)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  if (clusters == T) {
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  } else {
    sorted_labels <- paste(sort(levels(cluster_size$cluster),decreasing = T))
  }

  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  melt_mtx$tot_Ncells <-ave(melt_mtx$value, melt_mtx$variable, FUN=sum)
  melt_mtx$percent <- melt_mtx$value*100/melt_mtx$tot_Ncells
  colnames(melt_mtx)[2] <- "Dataset"

  # Plots per se
  if (value_reorder == T) {
    #Plot total proportion in Log scale.
    p1 <- ggplot(cluster_size, aes(x = value,y= reorder(cluster, value)))
    #Plot Barplot by ident
    p2 <- ggplot(melt_mtx,aes(x=reorder(cluster, percent),y=percent,fill=Dataset))
  } else {
    #Plot total proportion in Log scale.
    p1 <- ggplot(cluster_size, aes(x = value, y= cluster))
    #Plot Barplot by ident
    p2 <- ggplot(melt_mtx,aes(x=cluster,y=percent,fill=Dataset))
  }

  my_pal <- c(jdb_palette("corona"),jdb_palette("brewer_spectra"),jdb_palette("lawhoops"),
              jdb_palette("samba_color"),jdb_palette("wolfgang_basic"))

  p1 <- p1 + geom_bar(position="dodge", stat="identity",fill = "grey60") +
    theme_bw() + scale_x_log10() + xlab("Cells per ident, log10 scale") +
    {if(PlotNcells) geom_text(aes(y=cluster,x=1.1,label=value,hjust="bottom"), size=4)} + ylab("") +
    theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14))

  p2 <- p2 + geom_bar(position="dodge", stat="identity") + theme_bw() + coord_flip() +
    {if(PlotNcells) geom_text(aes(cluster,y = 0.1,label=value, group=Dataset, hjust =0),
                              position=position_dodge(width=1), size=4)} +
    scale_fill_manual(values = my_pal) +
    ylab("Percentage of cells per ident") + xlab("Idents") +
    theme(legend.position="top", axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14))

  #Sum plots
  p2 + p1 + plot_layout(widths = c(3,1))
}


