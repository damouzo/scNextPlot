

# Compile cellranger outputs ---------------------------------------------------
scPlot_CellRangerMultiSummary = function (base_path, secondary_path = "outs/",
                                          file_name="summary.csv", lib_list = NULL,
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
    raw_data <- read.csv(file = paste0(file_path, file_name),
                         stringsAsFactors = F)
    column_numbers <- grep(pattern = ",", x = raw_data[1,
    ])
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
  percent_cols <- grep("Percent|Valid|Fraction|Q30.bases|nuclear.read.pairs|mapped|with.TSO",
    colnames(full_data), value = TRUE)

  for (col in percent_cols) {
    full_data[[col]] <- full_data[[col]] * 100
  }

  full_data <- full_data %>%   rename_with(~ paste0(., ".percent"), all_of(percent_cols))



    # Add metadata extra
  if (!is.null(x = add_metadata)) {full_data <- merge(full_data, add_metadata, by="Sample.ID")}



  # Create beatufil table -----------------------------------------
  if (save_html) {
    df_all_plot <- full_data %>%
      mutate(Ncell_estim_cols = case_when(
        Estimated.number.of.cells <= 12000 & Estimated.number.of.cells >= 500 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(Feat_linkages_cols = case_when(Feature.linkages.detected >= 100 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_ValidBar_cols = case_when(ATAC.Valid.barcodes.percent >= 85 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Frac_genome_cols = case_when(ATAC.Fraction.of.genome.in.peaks.percent <= 75 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_TSS_cols = case_when(ATAC.TSS.enrichment.score >= 5 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Frac_over_peak_cols = case_when(ATAC.Fraction.of.high.quality.fragments.overlapping.peaks.percent >= 25 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Mean_reads = case_when(ATAC.Mean.raw.read.pairs.per.cell >= 5000 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Frac_cells_cols = case_when(ATAC.Fraction.of.high.quality.fragments.in.cells.percent >= 40 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Frac_peakscells_cols = case_when(ATAC.Fraction.of.transposition.events.in.peaks.in.cells.percent >= 25 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Median_frag_cols = case_when(ATAC.Median.high.quality.fragments.per.cell >= 100 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_Conf_Map_cols = case_when(ATAC.Confidently.mapped.read.pairs.percent >= 80 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(ATAC_NonNuclear_cols = case_when(ATAC.Non.nuclear.read.pairs.percent <= 10 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_ValidBar_cols = case_when(GEX.Valid.barcodes.percent >= 80 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_TSO_cols = case_when(GEX.Reads.with.TSO.percent <= 25 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Read_Map_Genome_cols = case_when(GEX.Reads.mapped.to.genome.percent >= 80 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Read_Intergenic_cols = case_when(GEX.Reads.mapped.confidently.to.intergenic.regions.percent <= 30 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Read_Transcriptome_cols = case_when(GEX.Reads.mapped.confidently.to.transcriptome.percent >= 50 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Read_Antisense_cols = case_when(GEX.Reads.mapped.antisense.to.gene.percent <= 30 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Mean_reads_cols = case_when(GEX.Mean.raw.reads.per.cell >= 5000 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Frac_Trans_cells_cols = case_when(GEX.Fraction.of.transcriptomic.reads.in.cells.percent >= 60 ~ "darkgreen", TRUE ~ "red"))%>%
      mutate(GEX_Median_UMI_cols = case_when(GEX.Median.UMI.counts.per.cell >= 100 ~ "darkgreen", TRUE ~ "red"))


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
                  "Estimated number of cells" = colDef(cell= data_bars(., text_position="above", text_color_ref="Ncell estim cols")),
                  "Feature linkages detected" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "Feat linkages cols")),
                  "ATAC Valid barcodes percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC ValidBar cols")),
                  "ATAC Fraction of genome in peaks percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Frac genome cols")),
                  "ATAC TSS enrichment score" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC TSS cols")),
                  "ATAC Fraction of high-quality fragments overlapping peaks percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Frac over peak cols")),
                  "ATAC Mean raw read pairs per cell" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Mean reads")),
                  "ATAC Fraction of high-quality fragments in cells percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Frac cells cols")),
                  "ATAC Fraction of transposition events in peaks in cells percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Frac peakscells cols")),
                  "ATAC Median high-quality fragments per cell" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Median frag cols")),
                  "ATAC Confidently mapped read pairs percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC Conf Map cols")),
                  "ATAC Non-nuclear read pairs percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "ATAC NonNuclear cols")),
                  "GEX Valid barcodes percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX ValidBar cols")),
                  "GEX Reads with TSO percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX TSO cols")),
                  "GEX Reads mapped to genome percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Read Map Genome cols")),
                  "GEX Reads mapped confidently to intergenic regions percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Read Intergenic cols")),
                  "GEX Reads mapped confidently to transcriptome percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Read Transcriptome cols")),
                  "GEX Reads mapped antisense to gene percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Read Antisense cols")),
                  "GEX Mean raw reads per cell" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Mean reads cols")),
                  "GEX Fraction of transcriptomic reads in cells percent" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Frac Trans cells cols")),
                  "GEX Median UMI counts per cell" = colDef(cell = data_bars(., text_position = "above", text_color_ref = "GEX Median UMI cols")),



                  # Remove for output
                  "Pipeline version" = colDef(show = FALSE),
                  "Ncell estim cols" = colDef(show = FALSE),
                  "Feat linkages cols" = colDef(show = FALSE),
                  "ATAC ValidBar cols" = colDef(show = FALSE),
                  "ATAC Frac genome cols" = colDef(show = FALSE),
                  "ATAC TSS cols" = colDef(show = FALSE),
                  "ATAC Frac over peak cols" = colDef(show = FALSE),
                  "ATAC Mean reads" = colDef(show = FALSE),
                  "ATAC Frac cells cols" = colDef(show = FALSE),
                  "ATAC Frac peakscells cols" = colDef(show = FALSE),
                  "ATAC Median frag cols" = colDef(show = FALSE),
                  "ATAC Conf Map cols" = colDef(show = FALSE),
                  "ATAC NonNuclear cols" = colDef(show = FALSE),
                  "GEX ValidBar cols" = colDef(show = FALSE),
                  "GEX TSO cols" = colDef(show = FALSE),
                  "GEX Read Map Genome cols" = colDef(show = FALSE),
                  "GEX Read Intergenic cols" = colDef(show = FALSE),
                  "GEX Read Transcriptome cols" = colDef(show = FALSE),
                  "GEX Read Antisense cols" = colDef(show = FALSE),
                  "GEX Mean reads cols" = colDef(show = FALSE),
                  "GEX Frac Trans cells cols" = colDef(show = FALSE),
                  "GEX Median UMI cols" = colDef(show = FALSE)

                )
      )
  }


  # Save files ---------------------------------------------
  if (save_html) {saveWidget(big_summary,title = "scPlotSummary",selfcontained= TRUE,libdir = "foldr2remove", file = paste0(output_path,"/scPlotSummary.html"))
    if (file.exists(paste0(output_path,"/foldr2remove"))) {unlink(paste0(output_path,"/foldr2remove"), recursive = TRUE)}}

  if (save_xlsx) {write.xlsx(full_data, file = paste0(output_path,"/scPlotSummary.xlsx"), sheetName = "Summary",row.names = F)}
  return(full_data)

}

