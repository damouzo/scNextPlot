
# Check Expression of Marker genes of SeuObj  ----------------------------------
scPlot_Heatmap <- function(SeuObj,
                           markers = NULL,
                           assay = DefaultAssay(SeuObj),
                           slot = "scaled.data",
                           need_scale = FALSE,
                           sort_var = NULL,
                           n = 8,
                           anno_var,
                           anno_colors,
                           hm_limit = NULL, #c(-2, 0, 2)
                           hm_colors = c("#4575b4","white","#d73027"),
                           cluster_rows = FALSE,
                           cluster_columns = TRUE,
                           show_column_names = FALSE,
                           show_row_names = TRUE,
                           row_names_side = "left",
                           show_column_dend = FALSE,
                           show_row_dend = FALSE,
                           row_font_size = 12
                           ) {

  # load required libraries-----------------------------------------------------
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(Signac)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(RColorBrewer)

  # Prepare Genes & Matrix -----------------------------------------------------
  # Gene set
  if (is.null(markers)) {
    markers <- SeuObj@assays[[assay]]@var.features
  }

  if (is.data.frame(markers)) {genes <- get_top_genes(SeuObj, markers, n)
  } else if (is.character(markers)) {genes <- markers
  } else {stop('Incorrect input of markers')}

  # Matrix

  mat <- as.matrix(GetAssayData(SeuObj, assay = assay, slot = slot)[markers, ])
  if (need_scale == T) {mat <- t(scale(t(mat)))}

  # Extract Annot Variables
  anno <- SeuObj@meta.data[, anno_var]


  #Work With Colors -----------------------------------------------------------
  #needed function
  are_colors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })
  }

  # Set Colors for variables side
  annos <- list()
  for (i in seq_along(1:length(anno_var))) {
    err_msg <- paste('Incorrect specification for annotation colors for', anno_var[i])
    value <- anno[[anno_var[i]]]

    if (is.numeric(value)) {
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != 'qual'])) {
        n <- brewer.pal.info[anno_colors[[i]],]['maxcolors'][[1]]
        pal <- brewer.pal(n = n, name = anno_colors[[i]])

        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)),
                              c(pal[2], pal[(n+1)/2], pal[n-1]))

      } else if (length(anno_colors[[i]]) == 3 & all(are_colors(anno_colors[[i]]))) {

        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), anno_colors[[i]])
      } else {
        stop(err_msg)
      }

      ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                              col = list(a = col_fun),
                              border = TRUE,
                              annotation_label = anno_var[i])
    } else {

      l <- levels(factor(anno[[anno_var[i]]]))
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
        col <- set_colors(anno_colors[[i]], length(l))

      } else if (length(anno_colors[[i]]) >= length(l) & all(are_colors(anno_colors[[i]]))) {
        col <- anno_colors[[i]]
      } else {
        stop(err_msg)
      }
      names(col) <- l
      col <- col[!is.na(names(col))]
      col <- list(a = col)

      ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                              col = col,
                              border = TRUE,
                              annotation_label = anno_var[i])
    }
    names(ha) <- anno_var[i]

    annos[[i]] <- ha
  }

  # Build annos sum up for all side colors
  annos <- do.call(c, annos)
  annos@gap <- rep(unit(1,"mm"), length(annos))
  ht_opt$message = FALSE

  # hm_limit
  if (is.null(hm_limit)) {
    hm_limit <- c(plyr::round_any(quantile(mat, c(0.1, 0.95))[[1]], 0.5, f=floor),
                  0,
                  plyr::round_any(quantile(mat, c(0.1, 0.95))[[2]], 0.5, f=ceiling))
  }

  col_fun = circlize::colorRamp2(hm_limit, hm_colors)

  # split and order
  if (!is.null(sort_var)) {
    factor_levels <-as.data.frame(as.matrix(table(anno[[sort_var]])))
    factor_levels <- row.names(factor_levels[order(-factor_levels$V1), , drop = FALSE])
    colum_split <- factor(anno[[sort_var]], levels = factor_levels)
    # HeatMap with order by var
    ht <- Heatmap(mat, column_split = colum_split,
                                column_title  = NULL,
                                cluster_rows = cluster_rows,
                                cluster_columns = cluster_columns,
                                heatmap_legend_param = list(direction = "horizontal",
                                                            legend_width = unit(6, "cm"),
                                                            title = "Expression"),
                                col = col_fun,
                                show_column_names = show_column_names,
                                show_row_names = show_row_names,
                                row_names_side = row_names_side,
                                show_column_dend = show_column_dend,
                                show_row_dend = show_row_dend,
                                row_names_gp = gpar(fontsize = row_font_size),
                                use_raster = T,raster_quality = 4,
                                top_annotation = annos)
  } else {
    # HeatMap withOUT order by var
    ht <- Heatmap(mat,
                  cluster_rows = cluster_rows,
                  cluster_columns = cluster_columns,
                  heatmap_legend_param = list(direction = "horizontal",
                                              legend_width = unit(6, "cm"),
                                              title = "Expression"),
                  col = col_fun,
                  show_column_names = show_column_names,
                  show_row_names = show_row_names,
                  row_names_side = row_names_side,
                  show_column_dend = show_column_dend,
                  row_names_gp = gpar(fontsize = row_font_size),
                  use_raster = T,raster_quality = 4,
                  top_annotation = annos)
  }

  draw(ht,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "right")

}











