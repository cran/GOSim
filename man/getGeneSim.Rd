\name{getGeneSim}
\alias{getGeneSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Compute functional similarity for genes  }
\description{
   Calculate the pairwise functional similarities for a list of genes
   using different strategies. 
}
\usage{
getGeneSim(genelist1, genelist2=NULL, similarity="funSimMax", similarityTerm="relevance", normalization=TRUE, method="sqrt", avg=(similarity=="OA"), verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist1}{ character vector of Entrez gene IDs }
  \item{genelist2}{ optional other character vector of Entrez gene IDs to compare against}
  \item{similarity}{ method to calculate the functional similarity between gene products }
  \item{similarityTerm}{ method to compute the similarity of GO terms }
  \item{normalization}{ normalize similarities yes/no}
  \item{method}{"sqrt": normalize sim(x,y) <- sim(x,y)/sqrt(sim(x,x)*sim(y,y)); "Lin": normalize sim(x,y) <- 2*sim(x,y)/(sim(x,x) + sim(y,y)); "Tanimoto": normalize sim(x,y) <- sim(x,y)/(sim(x,x) + sim(y,y) - sim(x,y)). NOTE: normalization does not have any effect, if similarity = "funSimMax", "funSimAvg" or similarity = "OA" and avg=TRUE}
  \item{avg}{ standardize the OA kernel by the maximum number of GO terms for both genes}
  \item{verbose}{ print out some information }
}
\details{
   The method to calculate the pairwise functional similarity between
   gene products can either be: 
   \describe{
\item{"max"}{the maximum similarity between any two GO terms}
\item{"mean"}{the average similarity between any two GO terms}
\item{funSimMax}{ the average of best matching GO term similarities. Take the maximum of the scores achieved by assignments of GO terms from gene 1 to gene 2 and vice versa. [2]}
\item{funSimAvg}{ the average of best matching GO term similarities. Take the average of the scores achieved by assignments of GO terms from gene 1 to gene 2 and vice versa. [2]}
\item{"OA"}{the optimal assignment (maximally weighted bipartite
  matching) of GO terms associated to the gene having fewer annotation
  to the GO terms of the other gene. [1]}
\item{"hausdorff"}{the translation of the Hausdorff distance between two sets: Let A and B be the two sets of GO terms associated to two genes. Then \eqn{sim(A,B) = min(min(max_{x \in A}(x,y)), min(max_{y\in B}(x,y)))}. [3]}
\item{"dot"}{the dot product between feature vectors describing the absence/presence of each GO term. The absence of each GO
 term is weighted by its information content. Depending on the type of later normalization one can arrive at the cosine similarity (method="sqrt") or at the Tanimoto coefficient (method="Tanimoto").[4]}
 }
}
\value{
  n x n similarity matrix (n = number of genes)
}
\references{ 
[1] H. Froehlich, N. Speer, C. Spieth, A. Zell, Kernel
  Based Functional Gene Grouping, Proc. Int. Joint Conf. on Neural
  Networks (IJCNN), 6886 - 6891, 2006.

[2] A. Schlicker, F. Domingues, J. Rahnenfuehrer, T. Lengauer, A new measure for functional similarity of gene products based on Gene Ontology, BMC Bioinformatics, 7, 302, 2006.

[3] A. del Pozo, F. Pazos, A. Valencia, Defining functional distances over Gene Ontology, BMC Bioinformatics, 9:50, 2008.

[4] M. Mistry, P Pavlidis, Gene Ontology term overlap as a measure of gene functional similarity, BMC Bioinformatics, 9:327, 2008.
} 
\author{ Holger Froehlich }
\note{ The result depends on the currently set ontology. 
}
\seealso{ \code{\link{getGeneSimPrototypes}}, \code{\link{getTermSim}}, \code{\link{setOntology}} }
\examples{
	# see evaluateClustering
}
\keyword{ file }% at least one, from doc/KEYWORDS
