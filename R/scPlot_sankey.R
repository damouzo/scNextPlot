
# Make a Sankey Plot From Seurat Meta.Data ################################################
scPlot_sankey <- function(SeuObj, steps, output_path=NULL) {

  # Load Libraries -------------------------------------------------
  set.seed(23)
  suppressMessages(library(Seurat))
  suppressMessages(library(dplyr))
  suppressMessages(library(networkD3))
  suppressMessages(library(htmlwidgets))

  # Settings -------------------------------------------------------
  if (is.null(output_path)) {output_path <- getwd()}

  # Set Data to plot -----------------------------------------------
  data <- SeuObj@meta.data[, steps]
  add_suffix <- function(x, suffix) {paste0(x, "_", suffix)}
  for (col in names(data)) {data[[col]] <- add_suffix(data[[col]], col)}


  # Create the df for the links in frecuency
  link_data <- data.frame(source = character(), target = character(), value = numeric(), stringsAsFactors = FALSE)

  for (i in 1:(length(steps)-1)) {
    freq_table <- as.data.frame(table(data[, c(steps[i], steps[i+1])]))
    names(freq_table) <- c("source", "target", "value")
    link_data <- rbind(link_data, freq_table)
  }


  # Data to sankey plot
  nodes <- unique(unlist(sapply(data, unique)))
  links <- data.frame(
    source = match(link_data$source, nodes) - 1,
    target = match(link_data$target, nodes) - 1,
    value = link_data$value
  )



  # Sankey Plot -----------------------------------------------------
  sankey_plot <- sankeyNetwork(
    Links = links,
    Nodes = data.frame(name = nodes),
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    fontSize = 13,
    sinksRight = FALSE
  )


  # Save HTML ---------------------------------------------------------
  saveWidget(sankey_plot, libdir = "foldr2remove", selfcontained = TRUE,
             file = paste0(output_path,"sankey_plot_", paste(colnames(data), collapse = "_"), ".html"))
  if (file.exists(paste0(output_path,"/foldr2remove"))) {unlink(paste0(output_path,"/foldr2remove"), recursive=TRUE)}
}


