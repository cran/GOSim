\name{getGeneSim}
\alias{getGeneSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Compute functional similarity for genes.  }
\description{
   Calculate the pairwise functional similarities for a list of genes
   using different strategies. 
}
\usage{
getGeneSim(genelist, similarity = "OA", similarityTerm = "JiangConrath", normalization = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{ character vector of Entrez gene IDs }
  \item{similarity}{ method to calculate the functional similarity between gene products }
  \item{similarityTerm}{ method to compute the similarity of GO terms }
  \item{normalization}{ normalize the similarities to [0,1] by transforming sim(x,y) <- sim(x,y)/sqrt(sim(x,x)*sim(y,y)) }
  \item{verbose}{ print out some information }
}
\details{
   The method to calculate the pairwise functional similarity between
   gene products can either be: 
\begin{ldescription}
\item{"max"}{the maximum similarity between any two GO terms}
\item{"OA"}{the optimal assignment (maximally weighted bipartite
  matching) of GO terms associated to the gene having fewer annotation
  to the GO terms of the other gene.}
\end{ldescription}
}
\value{
  n x n similarity matrix (n = number of genes)
}
\references{ H. Froehlich, N. Speer, C. Spieth, and A. Zell, Kernel
  Based Functional Gene Grouping, Proc. Int. Joint Conf. on Neural
  Networks (IJCNN), 6886 - 6891, 2006  } 
\author{ Holger Froehlich }
\note{ The result depends on the currently set ontology. 
}
\seealso{ \code{\link{getGeneSimPrototypes}}, \code{\link{getTermSim}}, \code{\link{setOntology}} }
\examples{
	# see evaluateClustering
}
\keyword{ file }% at least one, from doc/KEYWORDS
