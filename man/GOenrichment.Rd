\name{GOenrichment}
\alias{GOenrichment}
\alias{analyzeCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GO enrichment analysis }
\description{
  This function performs a GO enrichment analysis using topGO. It combines the two former functions "GOenrichment" and "analyzeCluster".
}
\usage{
GOenrichment(genesOfInterest, allgenes, cutoff=0.01, method=c("weightCount","elimCount"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
  \item{genesOfInterest}{character vector of Entrez gene IDs}
  \item{allgenes}{character vector of Entrez gene IDs}
  \item{cutoff}{significance cutoff for GO enrichment analysis}
  \item{method}{topGO method: either weight or elim}
}
\value{  
  \item{GOTerms}{list of significant GO terms and their description}
  \item{p.values}{vector of p-values for significant GO terms}
  \item{genes}{list of genes associated to each GO term}
}
\references{
Adrian Alexa, J\"org Rahnenf\"uhrer, Thomas Lengauer: Improved scoring of functional groups from gene expression data by decorrelating GO graph structure, Bioinformatics, 2006, 22(13):1600-1607}
\author{Holger Froehlich}
\seealso{\code{\link{evaluateClustering}} }
\examples{
	\dontrun{	
	setOntology("BP")
	gomap <- get("gomap",env=GOSimEnv)
	allgenes = sample(names(gomap), 1000) # suppose these are all genes
	allpvalues = runif(1000) # an these are their pvalues
	names(allpvalues) = allgenes
	if(require(topGO) & require(annotate))
		GOenrichment(allpvalues[allpvalues<0.05], allpvalues) # GO enrichment analysis
	}
}
\keyword{ file }% at least one, from doc/KEYWORDS
