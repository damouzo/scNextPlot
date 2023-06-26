

# Compile cellranger outputs ---------------------------------------------------
scPlot_CellRangerCiteSummary = function (base_path, secondary_path = "outs/",
                                          file_name="metrics_summary.csv", lib_list = NULL,
                                          new_lib_names = NULL, add_metadata= NULL,
                                          output_path=NULL, save_html=T, save_xlsx=T) {


  # Load Libraries -------------------------------------------------
  suppressMessages(library(Seurat))
  suppressMessages(library(pbapply))
  suppressMessages(library(dplyr))
  suppressMessages(library(reactable))
  suppressMessages(library(reactablefmtr))
  suppressMessages(library(htmlwidgets))
  suppressMessages(library(cli))
  suppressMessages(library(xlsx))



  # Take Summary of every sample ------------------------------------
  # Based on scCustomize
  if (dir.exists(paths = base_path) == FALSE) {cli_abort(message = "Directory base_path does not exist.")}
  if (is.null(x = lib_list)) {lib_list <- list.dirs(path = base_path, full.names = F, recursive = F)}
  if (is.null(output_path)) {output_path <- getwd()}
  if (is.null(x = secondary_path)) {secondary_path <- ""}

  for (i in 1:length(x = lib_list)) {
    full_directory_path <- file.path(base_path, lib_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {cli_abort(message = "Directory under basepath not exist")}}

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    }
    else {
      file_path <- file.path(base_path, lib_list[x], secondary_path)
    }
    raw_data <- read.csv(file= paste0(file_path, file_name),stringsAsFactors= F)
    column_numbers <- grep(pattern = ",", x = raw_data[1,])
    raw_data[, c(column_numbers)] <- lapply(raw_data[, c(column_numbers)],
                                            function(x) {as.numeric(gsub(",", "", x))})
    return(raw_data)
  })

  if (is.null(x = new_lib_names)) {
    names(raw_data_list) <- lib_list
  }  else {names(raw_data_list) <- new_lib_names}

  full_data <- bind_rows(raw_data_list, .id = "sample_id")
  colnames(full_data) <- gsub(pattern = "\\.", replacement = "_",
                              x = colnames(x = full_data))
  colnames(full_data) <- gsub("_", ".", colnames(full_data))
  full_data <- full_data[, -1] #remove doble ID column

  # Add percent format
  percent_cols <- colnames(full_data)[apply(full_data, 2, function(x) any(grepl("%", x)))]
  full_data <- full_data %>%   rename_with(~ paste0(., ".percent"), all_of(percent_cols))

  for (col in colnames(full_data)) {
    if (is.character(full_data[[col]])) {
      full_data[[col]] <- as.numeric(gsub("%", "", full_data[[col]]))
    }
  }

  # Add Sample ID
  full_data$Sample.ID <- lib_list
  full_data <- full_data[, c("Sample.ID", setdiff(colnames(full_data), "Sample.ID"))]

  # Add metadata extra
  if (!is.null(x = add_metadata)) {full_data <- merge(full_data, add_metadata, by="Sample.ID")}



  # Create beatufil table -----------------------------------------
  if (save_html) {
    df_all_plot <- full_data %>%
      mutate(Ncell_estim_cols = case_when(
        Estimated.Number.of.Cells <= 12000 & Estimated.Number.of.Cells >= 500 ~ "darkgreen", TRUE ~ "red"))#%>%
      # mutate(Feat_linkages_cols = case_when(Feature.linkages.detected >= 100 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_ValidBar_cols = case_when(ATAC.Valid.barcodes.percent >= 85 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Frac_genome_cols = case_when(ATAC.Fraction.of.genome.in.peaks.percent <= 75 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_TSS_cols = case_when(ATAC.TSS.enrichment.score >= 5 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Frac_over_peak_cols = case_when(ATAC.Fraction.of.high.quality.fragments.overlapping.peaks.percent >= 25 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Mean_reads = case_when(ATAC.Mean.raw.read.pairs.per.cell >= 5000 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Frac_cells_cols = case_when(ATAC.Fraction.of.high.quality.fragments.in.cells.percent >= 40 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Frac_peakscells_cols = case_when(ATAC.Fraction.of.transposition.events.in.peaks.in.cells.percent >= 25 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Median_frag_cols = case_when(ATAC.Median.high.quality.fragments.per.cell >= 100 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_Conf_Map_cols = case_when(ATAC.Confidently.mapped.read.pairs.percent >= 80 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(ATAC_NonNuclear_cols = case_when(ATAC.Non.nuclear.read.pairs.percent <= 10 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_ValidBar_cols = case_when(GEX.Valid.barcodes.percent >= 80 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_TSO_cols = case_when(GEX.Reads.with.TSO.percent <= 25 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Read_Map_Genome_cols = case_when(GEX.Reads.mapped.to.genome.percent >= 80 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Read_Intergenic_cols = case_when(GEX.Reads.mapped.confidently.to.intergenic.regions.percent <= 30 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Read_Transcriptome_cols = case_when(GEX.Reads.mapped.confidently.to.transcriptome.percent >= 50 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Read_Antisense_cols = case_when(GEX.Reads.mapped.antisense.to.gene.percent <= 30 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Mean_reads_cols = case_when(GEX.Mean.raw.reads.per.cell >= 5000 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Frac_Trans_cells_cols = case_when(GEX.Fraction.of.transcriptomic.reads.in.cells.percent >= 60 ~ "darkgreen", TRUE ~ "red"))%>%
      # mutate(GEX_Median_UMI_cols = case_when(GEX.Median.UMI.counts.per.cell >= 100 ~ "darkgreen", TRUE ~ "red"))


    colnames(df_all_plot) <- gsub("\\.", " ", colnames(df_all_plot))
    colnames(df_all_plot) <- gsub("_", " ", colnames(df_all_plot))


    big_summary <- df_all_plot %>%
      reactable(.,
                theme = fivethirtyeight(centered = TRUE, header_font_size = 11),     # or nytimes(),fivethirtyeight
                defaultSorted = "Sample ID",
                searchable = F,
                resizable = TRUE,
                defaultColDef = colDef(cell = data_bars(., text_size = 13, text_position = "above"),searchable = TRUE,
                                       minWidth = 160, align = "center", vAlign = "center",headerVAlign="center"),

                columns = list(
                  # Section div.
                  "Sample ID" = colDef(sticky = "left", style=list(background = "#bdbdbd",borderRight = "1px solid #bdbdbd"),
                                       headerStyle=list(background="#bdbdbd",borderRight="1px solid #bdbdbd",color = "black"), minWidth = 100),

                  # Fix Lengh
                  "Genome" = colDef(minWidth = 90),

                  # Add color levels
                  "Estimated Number of Cells" = colDef(cell= data_bars(., text_position="above", text_color_ref="Ncell estim cols")),


                  # Remove for output
                  "Pipeline version" = colDef(show = FALSE),
                  "Ncell estim cols" = colDef(show = FALSE)
                )
      )
  }


  # Save files ---------------------------------------------
  if (save_html) {saveWidget(big_summary,title = "scPlotSummary",selfcontained= TRUE,libdir = "foldr2remove", file = paste0(output_path,"/scPlotSummary.html"))
    if (file.exists(paste0(output_path,"/foldr2remove"))) {unlink(paste0(output_path,"/foldr2remove"), recursive = TRUE)}}

  if (save_xlsx) {write.xlsx(full_data, file = paste0(output_path,"/scPlotSummary.xlsx"), sheetName = "Summary",row.names = F)}
  return(full_data)

}

