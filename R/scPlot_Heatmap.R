
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
                           row_col_var_df = NULL,
                           row_col_var_plot = NULL,
                           hm_limit = NULL, #c(-2, 0, 2)
                           hm_colors = c("#4575b4","white","#d73027"),
                           hm_scale_name = "Expression",
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

  # Set Colors for variables Top  ----------------------------------------------
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

  # Build annos sum up for all side colors -------------------------------------
  annos <- do.call(c, annos)
  annos@gap <- rep(unit(1,"mm"), length(annos))
  ht_opt$message = FALSE

  # hm_limit ------------------------------------------------------------------
  if (is.null(hm_limit)) {
    hm_limit <- c(plyr::round_any(quantile(mat, c(0.1, 0.95))[[1]], 0.5, f=floor),
                  0,
                  plyr::round_any(quantile(mat, c(0.1, 0.95))[[2]], 0.5, f=ceiling))
  }

  col_fun = circlize::colorRamp2(hm_limit, hm_colors)

  # Set row annot Colors for variables ----------------------------------------
  if(!is.null(row_col_var_df)){
    # obtain only column for plots
    row_col_var_df <- row_col_var_df %>% dplyr::select(all_of(row_col_var_plot))

    # Set colors function
    generate_color_palette <- function(df) {
      color_palettes <- list()

      # Our paletes for row_cols
      numeric_gradients <- list(c("white", "orange"), c("white", "purple"),
                    c("white", "blue"), c("white", "green"), c("white", "red"))
      categorical_palettes <- list("Set3","Accent","Paired", "Dark2", "Set2")
      numeric_gradient_index <- 1
      categorical_palette_index <- 1


      for (column_name in colnames(df)) {
        unique_items <- unique(df[[column_name]])
        num_unique_items <- length(unique_items)

        if (is.numeric(df[[column_name]])) {# Paleta continua para variables numéricas
          col_func <- colorRamp2(c(min(df[[column_name]], na.rm = TRUE),
                                   max(df[[column_name]], na.rm = TRUE)),
                                 numeric_gradients[[numeric_gradient_index]])
          color_palettes[[column_name]] <- col_func
          numeric_gradient_index <- numeric_gradient_index %% length(numeric_gradients) + 1

        } else { # Paleta discreta para variables categóricas
          pal <- brewer.pal(min(num_unique_items, 9), categorical_palettes[[categorical_palette_index]])
          color_palettes[[column_name]] <- setNames(pal[1:num_unique_items], unique_items)
          categorical_palette_index <- categorical_palette_index %% length(categorical_palettes) + 1

        }
      }

      return(color_palettes)
    }
    row_col_pal <- generate_color_palette(row_col_var_df)

    left_annotation <- rowAnnotation(df = row_col_var_df, col = row_col_pal)
  } else {
    left_annotation <- NULL
  }


  # HEATMAPS ###################################################################
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
                                                            title = hm_scale_name),
                                col = col_fun,
                                show_column_names = show_column_names,
                                left_annotation=left_annotation,
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
                                              title = hm_scale_name),
                  col = col_fun,
                  show_column_names = show_column_names,
                  show_row_names = show_row_names,
                  left_annotation=left_annotation,
                  row_names_side = row_names_side,
                  show_column_dend = show_column_dend,
                  row_names_gp = gpar(fontsize = row_font_size),
                  use_raster = T,raster_quality = 4,
                  top_annotation = annos)
  }

  draw(ht, heatmap_legend_side="bottom", annotation_legend_side="right",legend_grouping="original")

}




