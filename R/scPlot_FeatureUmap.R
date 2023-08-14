

# Check Expression of Marker genes of SeuObj  ----------------------------------
scPlot_FeatureUmap <- function(obj = countData, features = "nGene", feature.type = "meta", umap = "umap",
                              ncols = ceiling(sqrt(length(features))), pt.size = 1, same.scale = FALSE, title = "GEX",
                              lower = NULL, upper = NULL, na.color = "gray") {

  # load required libraries-----------------------------------------------------
  library(reshape2)
  library(plotly)
  library(tidyr)
  library(gridExtra)

  colors <- c("#191970", "#121285", "#0C0C9A", "#0707B0", "#0101C5", "#0014CF", "#0033D3", "#0053D8", "#0072DD",
              "#0092E1", "#00B2E6", "#00D1EB", "#23E8CD", "#7AF17B", "#D2FA29", "#FFEB00", "#FFC300", "#FF9B00",
              "#FF8400", "#FF7800", "#FF6B00", "#FF5F00", "#FF5300", "#FF4700", "#F73B00", "#EF2E00", "#E62300",
              "#DD1700", "#D50B00", "#CD0000")
  ncols=ceiling(sqrt(length(features)))

  if (feature.type != "meta" & feature.type != "gene") {
    stop("feature type must be 'meta' or 'gene'")
  }

  # Print Meta Feature ---------------------------------------------------------
  if (feature.type == "meta") {
    feature.info <- as.matrix(obj@meta.data[, features])
    if (length(features) == 1) {
      colnames(feature.info) <- features
    }

    tmp.df <- data.frame(feature.info, obj@reductions[[umap]]@cell.embeddings)
    plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key = TRUE)

  } else if (feature.type == "gene") {
    featuresNotFound <- features[!features %in% rownames(obj)]
    if (length(featuresNotFound) != 0) {print(paste0("Gene(s) ",featuresNotFound,
                                                     " not found in Seurat Object" ))}
    features <- features[features %in% rownames(obj)]
    feature.info <- as.matrix(GetAssayData(object = obj, slot = "data")[features,])
    if (length(features) == 1) {
      colnames(feature.info) <- features
      feature.info <- t(feature.info)
    }

    tmp.df <- data.frame(t(feature.info), obj@reductions[[umap]]@cell.embeddings)
    plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key = TRUE)
  }

  if (!same.scale) {
    p_list <- list()

    for (i in 1:length(features)) {
      feature_data <- plot.df[plot.df$name == features[i], ]
      feature_lower <- min(feature_data$val)
      feature_upper <- max(feature_data$val)

      p <- ggplot(feature_data, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = val, alpha = val), shape = 19, size = pt.size) +
        theme(aspect.ratio = 1) +
        scale_color_gradientn(colors = colors, limits = c(feature_lower, feature_upper),
                              na.value = na.color) +
        scale_alpha_continuous(range = c(0.6, 1), guide = "none") +
        labs(color = features[i], title = title) +
        theme(
          aspect.ratio = 1,
          text = element_text(size = 10),
          axis.text = element_text(size = 6),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
          strip.text = element_text(size = 10),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(linewidth = 0.5)
        )

      p_list[[i]] <- p
    }

    #combined_plot <- do.call(grid.arrange, c(p_list, ncol = ncols))
    combined_plot <- cowplot::plot_grid(plotlist = p_list, ncol = ncols)
    return(combined_plot)

  } else {
    if (!length(lower) & !length(upper)) {
      lower <- min(plot.df$val)
      upper <- max(plot.df$val)
    }

    p <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = val,alpha = val), shape = 19, size = pt.size) +
      theme(aspect.ratio = 1) +
      scale_color_gradientn(colors = colors, limits = c(lower, upper), na.value = na.color) +
      scale_alpha_continuous(range = c(0.6, 1), guide = "none") +
      labs(color = title, title = title) +
      facet_wrap(~name, ncol = ncols) +
      theme(
        aspect.ratio = 1,
        text = element_text(size = 10),
        axis.text = element_text(size = 6),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
        strip.text = element_text(size = 10),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(linewidth = 0.5)
      )

    return(p)
  }
}



