\name{getGeneSim}
\alias{getGeneSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Compute functional similarity for genes.  }
\description{
   Calculate the pairwise functional similarities for a list of genes
   using different strategies. 
}
\usage{
getGeneSim(genelist, similarity = "funSimMax", similarityTerm = "Lin", normalization = TRUE, method = "sqrt", avg = (similarity=="OA"), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{ character vector of Entrez gene IDs }
  \item{similarity}{ method to calculate the functional similarity between gene products }
  \item{similarityTerm}{ method to compute the similarity of GO terms }
  \item{normalization}{ normalize similarities yes/no}
  \item{method}{"sqrt": normalize sim(x,y) <- sim(x,y)/sqrt(sim(x,x)*sim(y,y)); "Lin": normalize sim(x,y) <- 2*sim(x,y)/(sim(x,x) + sim(y,y)). NOTE: normalization does not have any effect, if similarity = "funSimMax", "funSimAvg" or similarity = "OA" and avg=TRUE}
  \item{avg}{ standardize the OA kernel by the maximum number of GO terms for both genes}
  \item{verbose}{ print out some information }
}
\details{
   The method to calculate the pairwise functional similarity between
   gene products can either be: 
\begin{ldescription}
\item{"max"}{the maximum similarity between any two GO terms}
\item{"mean"}{the average similarity between any two GO terms}
\item{funSimMax}{ the average of best matching GO term similarities. Take the maximum of the scores achieved by assignments of GO terms from gene 1 to gene 2 and vice versa.}
\item{funSimAvg}{ the average of best matching GO term similarities. Take the average of the scores achieved by assignments of GO terms from gene 1 to gene 2 and vice versa.}
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
  Networks (IJCNN), 6886 - 6891, 2006.

A. Schlicker, F. Domingues, J. Rahnenfuehrer and Thomas Lengauer, A new measure for functional similarity of gene products based on Gene Ontology, BMC Bioinformatics, 7, 302, 2006.
} 
\author{ Holger Froehlich }
\note{ The result depends on the currently set ontology. 
}
\seealso{ \code{\link{getGeneSimPrototypes}}, \code{\link{getTermSim}}, \code{\link{setOntology}} }
\examples{
	# see evaluateClustering
}
\keyword{ file }% at least one, from doc/KEYWORDS
