\name{getGeneFeaturesPrototypes}
\alias{getGeneFeaturesPrototypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get feature vector representation of genes. }
\description{
  Computes the feature vectors for list of genes: Each gene is
  represented by its similarities to predefined prototype genes. 
}
\usage{
getGeneFeaturesPrototypes(genelist, prototypes = NULL,
                          similarity = "max", similarityTerm = "JiangConrath",
                          pca = TRUE, normalization = TRUE, verbose = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{ character vector of Entrez gene IDs }
  \item{prototypes}{ character vector of Entrez gene IDs representing the prototypes }
  \item{similarity}{ method to calculate the similarity to prototypes }
  \item{similarityTerm}{ method to compute the GO term similarity }
  \item{pca}{ perform PCA on feature vectors to reduce dimensionality }
  \item{normalization}{ scale the feature vectors to norm 1}
  \item{verbose}{ print out additional information }
}
\details{
If no prototypes are passed, the method calls the
\code{\link{selectPrototypes}} function with no arguments. Hence, the
prototypes in this case are the 250 genes with most known annotations.

The PCA postprocessing determines the principal components that can
explain at least 95\% of the total variance in the feature space. 

The method to calculate the functional similarity of a gene to a certain prototype can either be
\item{"max"}{the maximum similarity between any two GO terms}
\item{"OA"}{the optimal assignment (maximally weighted bipartite
  matching) of GO terms associated to the gene having fewer annotation
  to the GO terms of the other gene.}
}
\value{
List with items
  \item{"features"}{feature vectors for each gene: n x d data matrix}
  \item{"prototypes"}{prototypes (= prinicipal components, if PCA has been performed)}  
}
\references{ [1] H. Froehlich, N. Speer, C. Spieth, and A. Zell, Kernel Based Functional Gene Grouping, Proc. Int. Joint Conf. on Neural Networks (IJCNN), 6886 - 6891, 2006 \newline
[2] N. Speer, H. Froehlich, A. Zell, Functional Grouping of Genes Using Spectral Clustering and Gene Ontology, Proc. Int. Joint Conf. on Neural Networks (IJCNN), pp. 298 - 303, 2005 }
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

