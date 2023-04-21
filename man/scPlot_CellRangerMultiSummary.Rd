\name{scPlot_CellRangerMultiSummary}
\alias{scPlot_CellRangerMultiSummary}
\title{scPlot_CellRangerMultiSummary}
\description{
This is a function that compiles the summary outputs from Cell Ranger and generates a user-friendly table with estimated cell counts. It also allows for the option to save the table as an HTML file or an Excel file.
}
\usage{
scPlot_CellRangerMultiSummary(base_path = path_counts, output_path = path_plot)
}
\arguments{
  \item{base_path}{Path to folder over out/ folder of Cellranger count function.}
  \item{output_path}{Path where save the output file.}
  \item{secondary_path}{If not default folder name out/ set new folder name. [Default:NULL]}
  \item{file_name}{Name of the file with stored summary. [Default:"summary.csv"]}
  \item{lib_list}{Library names to only check that library from. [Default: All Samples] }
  \item{new_lib_names}{Vector with same size of lib_list with library titles to rename. [Default:NULL]}
  \item{save_html}{Bolean for save or not HTML reactable file. [Default:T]}
  \item{save_xlsx}{Bolean for save or not XLSX file. [Default:T]}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Daniel Mouzo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  # After running CellRanger
path_counts <- "path/to/folder/cellranger/Count/output/"
path_plot <- "path/to/save/output/summary"
#secondary_path by default is "outs/"

scPlot_CellRangerMultiSummary(base_path = path_counts, output_path = path_plot,
                              save_html = TRUE, save_xlsx = TRUE)


{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
