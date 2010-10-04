\name{getGeneSimPrototypes}
\alias{getGeneSimPrototypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute functional similarity of genes with respect to a feature vector representation. }
\description{
  Computes the pairwise functional similarities for a list of genes: Each
  gene is represented by a feature vector containing the gene's
  similarities to predefined prototype genes. 
}
\usage{
getGeneSimPrototypes(genelist, prototypes = NULL, similarity = "max",
                     similarityTerm = "JiangConrath", pca = TRUE,
                     normalization = TRUE, verbose = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{ character vector of Entrez gene IDs }
  \item{prototypes}{ character vector of Entrez gene IDs representing the prototypes }
  \item{similarity}{ method to calculate the similarity to prototypes }
  \item{similarityTerm}{ method to compute the GO term similarity }
  \item{pca}{ perform PCA on feature vectors to reduce dimensionality }
  \item{normalization}{ normalize similarities to [0,1]: sim(x,y) <- 0.5*(sim(x,y)/sqrt(sim(x,x)*sim(y,y)) + 1) } 
  \item{verbose}{print additional information}
}
\details{
  The method calls \code{\link{getGeneFeaturesPrototypes}} to calculate
  the feature vectors. The functional similarity between two genes is
  essentially given by the dot product between their feature vectors. 
}
\value{
List with items
  \item{"similarity"}{n x n similarity matrix (n = number of genes)}
  \item{"prototypes"}{prototypes (= prinicipal components, if PCA has been performed)}
  \item{"features"}{feature vectors for each gene: n x d data matrix}
}
\references{ 
[1] H. Froehlich, N. Speer, C. Spieth, and A. Zell, Kernel Based Functional Gene Grouping, Proc. Int. Joint Conf. on Neural Networks (IJCNN), 6886 - 6891, 2006

[2] N. Speer, H. Froehlich, A. Zell, Functional Grouping of Genes Using Spectral Clustering and Gene Ontology, Proc. Int. Joint Conf. on Neural Networks (IJCNN), pp. 298 - 303, 2005 }
\author{ Holger Froehlich }
\note{ The result depends on the currently set ontology ("BP","MF","CC").
}
\seealso{ \code{\link{getGeneFeaturesPrototypes}},
  \code{\link{selectPrototypes}}, \code{\link{getGeneSim}},
  \code{\link{getTermSim}}, \code{\link{setOntology}}  }
\examples{
 \dontrun{ may take some time ...
 proto=selectPrototypes(n=50) # --> returns a character vector of 50 genes with the highest number of annotations
 getGeneSimPrototypes(c("207","208","360"),prototypes=proto, similarityTerm="Resnik")
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
