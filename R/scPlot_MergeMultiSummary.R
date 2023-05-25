


# Join Summaries of cellranger outputs -----------------------------------------
scPlot_MergeMultiSummary = function (xlsx2merge,output_path=NULL,
                                     save_html=T, save_xlsx=T) {


  # Load Libraries -------------------------------------------------
  suppressMessages(library(Seurat))
  suppressMessages(library(pbapply))
  suppressMessages(library(dplyr))
  suppressMessages(library(reactable))
  suppressMessages(library(reactablefmtr))
  suppressMessages(library(htmlwidgets))
  suppressMessages(library(cli))
  suppressMessages(library(xlsx))

  # Load Summary of every sample ----------------------------------
  AllSumaries <- list()
  for (i in 1:length(xlsx2merge)) {
    file <- xlsx2merge[i]
    data <- read.xlsx(file, sheetIndex = 1)

    AllSumaries[[i]] <- data
  }

  # Merge ---------------------------------------------------------
  full_data <- bind_rows(AllSumaries)



  # Create beatufil table -----------------------------------------
  if (save_html) {
    df_all_plot <- full_data %>%
      mutate(Ncell_estim_cols = case_when(Estimated.number.of.cells <= 11000 ~ "darkgreen", TRUE ~ "red"))

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

                  # Remove for output
                  "Pipeline version" = colDef(show = FALSE),
                  "Ncell estim cols" = colDef(show = FALSE)
                )
      )
  }


  # Save files ---------------------------------------------
  if (save_html) {saveWidget(big_summary,title = "scPlotSummaryMerged",selfcontained= TRUE,libdir = "foldr2remove", file = paste0(output_path,"/scPlotSummaryMerged.html"))
    if (file.exists(paste0(output_path,"/foldr2remove"))) {unlink(paste0(output_path,"/foldr2remove"), recursive = TRUE)}}

  if (save_xlsx) {write.xlsx(full_data, file = paste0(output_path,"/scPlotSummaryMerged.xlsx"), sheetName = "Summary",row.names = F)}
  return(full_data)

}
