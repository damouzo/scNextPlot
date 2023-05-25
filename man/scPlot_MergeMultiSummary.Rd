\name{scPlot_MergeMultiSummary}
\alias{scPlot_MergeMultiSummary}
\title{scPlot_MergeMultiSummary - Function to merge multiple summaries of cellranger outputs}
\description{
This function merges multiple summaries of cellranger outputs stored in XLSX files. It combines the rows of each summary, with missing columns filled with NA values.
}
\usage{
scPlot_MergeMultiSummary(xlsx2merge, output_path = NULL, save_html = TRUE, save_xlsx = TRUE)
}
\arguments{
  \item{xlsx2merge}{
    A character vector specifying the paths to the XLSX files containing the summaries to be merged.
  }

  \item{output_path}{
    Optional. The path to the output directory where the merged summary files will be saved. If not provided, the files will be saved in the current working directory.
  }

  \item{save_html}{
    Optional. A logical value indicating whether to save the merged summary as an HTML file. Default is TRUE.
  }

  \item{save_xlsx}{
    Optional. A logical value indicating whether to save the merged summary as an XLSX file. Default is TRUE.
  }
}
\details{
The function reads each XLSX file from the paths specified in \code{xlsx2merge}, loads the summaries into memory, and then merges them into a single data frame. The resulting merged summary can be saved as an HTML file and/or an XLSX file.
}
\value{
The function returns the merged summary as a data frame. Additionally, it saves the merged summary as an HTML file and/or an XLSX file if specified.
}
\author{
Daniel Mouzo
}
\note{
This function requires the following packages to be installed: \code{Seurat}, \code{pbapply}, \code{dplyr}, \code{reactable}, \code{reactablefmtr}, \code{htmlwidgets}, \code{cli}, and \code{xlsx}.
}

\seealso{
\code{\link{read.xlsx}}, \code{\link{bind_rows}}, \code{\link{saveWidget}}, \code{\link{write.xlsx}}
}

\examples{
## Merge summaries of two XLSX files and save the merged summary as an HTML and XLSX file
xlsx2merge <- c("path/to/first.xlsx", "path/to/second.xlsx")
scPlot_MergeMultiSummary(xlsx2merge, output_path = "path/to/output")
}
