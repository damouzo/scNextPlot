\name{scPlot_allres}
\alias{scPlot_allres}
\title{scPlot_allres}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scPlot_allres(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Seurat Object
}
 \item{res}{
Select the resolutions to test [Default: c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)].
}
 \item{clustree}{
True for get the clustertree with all the tested resolutions. [Default:T]
}
 \item{UMAPs}{
True for get the individual UMAP of each resolution. [Default:T]
}
 \item{clustreeGenes}{
Vector of genes to check the expression in the Clustertree. [Default:NULL]
}
 \item{reduction}{
Seleccion of the reduction for the DimPlot. [Default:"umap"]
}
 \item{Graph}{
Custom name of the graph made by FindNeighbors. [Default:NULL]
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
  # Step Before
Seurat_Object <- FindNeighbors(Seurat_Object, dims = 1:30, verbose = F)

  # Plot selected resolution test
pdf("AllRes.pdf", width=10, height=10)
scPlot_allres(integrated_sct, res=c(0,0.5,0.6,0.7,1), reduction="umap",
              clustree = T,UMAPs=T, clustreeGenes = c("LYZ"))
dev.off()

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
