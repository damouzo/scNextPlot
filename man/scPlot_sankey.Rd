\name{scPlot_sankey}
\alias{scPlot_sankey}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{scPlot_sankey
%%  ~~function to do ... ~~
}
\description{
Create a Sankey plot to check meta.data label, for example for different annotations.
}
\usage{
scPlot_sankey(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{SeuObj}{
Seurat Object
}
 \item{steps}{
names of the 2 columns of meta.data that want to check (e.g. c("Azimuth"", "SingleR""))
}
 \item{output_path}{
Paths where save the html [Default: getwd()]
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
scPlot_sankey(SeuObj, nodes = c("SingleR.fine", "LabelTransfer"),
                output_path="/path/for/save/plot/")
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
