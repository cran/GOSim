\name{getGeneFeatures}
\alias{getGeneFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get simple feature vector representation of genes}
\description{
  Computes feature vectors for list of genes: Each gene is
  represented by a vector describing the presence/absence of all GO terms. The absence of each GO term is additionally weighted by its information content. 
}
\usage{
getGeneFeatures(genelist, pca=FALSE, normalization=FALSE, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{ character vector of Entrez gene IDs }
  \item{pca}{ perform PCA on feature vectors to reduce dimensionality }
  \item{normalization}{ scale the feature vectors to norm 1}
  \item{verbose}{ print out additional information }
}
\details{
The PCA postprocessing determines the principal components that can
explain at least 95\% of the total variance in the feature space. 
}
\value{
matrix with rows being genes and columns being GO terms.
}
\references{ 
M. Mistry, P Pavlidis, Gene Ontology term overlap as a measure of gene functional similarity, BMC Bioinformatics, 9:327, 2008.
}
\author{ Holger Froehlich }
\note{ The result depends on the currently set ontology ("BP","MF","CC").
}
\seealso{ \code{\link{getGeneSimPrototypes}},
  \code{\link{selectPrototypes}}, \code{\link{getGeneSim}},
  \code{\link{getTermSim}}, \code{\link{setOntology}} } 
\examples{
	# see selectPrototypes
}
\keyword{ file }% at least one, from doc/KEYWORDS

