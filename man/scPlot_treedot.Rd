\name{scPlot_treedot}
\alias{scPlot_treedot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
scPlot_treedot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scPlot_treedot(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object with multiple GSEA data from pairwise_termsim()
}
  \item{top_paths}{
# of how many paths of each compression take into account. [Default:5]
}
  \item{clust_num}{
# of clusters into which you want to separate the dendrogram. [Default:3]
}
  \item{ORA_type}{
ORA function selection ("GO", "KEGG", "DO"). [Default:"GO"]
}
  \item{ORA_ont}{
What ontology do you want to use? [Default:"BP"]
}
  \item{ORA_minGSSize}{
The min. size of genes annotated by Ontology term for testing. [Default:10]
}
  \item{ORA_maxGGSSize}{
The max. size of genes annotated by Ontology term for testing. [Default:500]
}
  \item{ORA_p.adj}{
P-value adjust cutoff.[Default:1]
}
  \item{ORA_p.adj_Meth}{
Method for correction by multiple tests for ORA [Default:"BH"]
}
  \item{ORA_KEGG_Org}{
Organism selection for the ORA performed with KEGG [Default:"hsa"]
}
  \item{ORA_GO_OrgDb}{
Data base selection for the ORA performed with GO [Default:"org.Hs.eg.db"]
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
compare_KEGG2 <- compareCluster(geneClusters=compare_KEGG, fun="gseKEGG",
                                organism="human",keyType="kegg",seed=23,
                                minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)
compare_KEGG2 <- pairwise_termsim(compare_KEGG2)


plot_treedot(comp_pair_term, top_paths=5, clust_num=3, ORA_ont="BP",
             ORA_minGSSize=10, ORA_maxGGSSize=500, ORA_p.adj=1)

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
