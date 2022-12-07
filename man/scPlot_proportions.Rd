\name{scPlot_proportions}
\alias{scPlot_proportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
scPlot_proportions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scPlot_proportions(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Seurat Object
}
  \item{by_ident}{
Select the categorical column [of SeuObj@meta.data] that interests you.
}
  \item{clusters}{
If the ident are cluster -> T, if it is labels -> F. [Default:F]
}
  \item{PlotNcells}{
Plot the number of cells within that category. [Default:T]
}
  \item{value_reorder}{
True to sort categories by number or alphabetically. [Default:T]
}
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
  # Plot Clusters Proportions
scPlot_proportions(seurat_obj, seurat_clusters, clusters=T,
                   PlotNcells=T, value_reorder=T)

  # Plot Label Proportions
scPlot_proportions(seurat_obj, predicted.celltype.l2, clusters=F,
                   PlotNcells=T, value_reorder=T)
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
