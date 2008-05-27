\name{selectPrototypes}
\alias{selectPrototypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Heuristic selection of prototypes and dimensionality reduction of feature vectors. }
\description{
\item Heuristic selection of prototypes 
\item Dimensionality reduction of feature vectors
}
\usage{
selectPrototypes(n = 250, method = "frequency", data = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of prototypes or maximum number of clusters }
  \item{method}{ method to select prototypes or to perform subset selection }
  \item{data}{ data matrix (l x d) of feature vectors (l = number of genes) }
  \item{verbose}{ print out information }
}
\details{
The following heuristics to perform automatic selection of prototypes are implemented:
\begin{ldescription}
\item{"frequency"}{ select n genes with highest number of GO annotations in the currently selected ontology}
\item{"random"}{ select n genes uniform randomly over all genes with annotations in the currently selected ontology}
\end{ldescription}

To perfom dimensionality reduction implemented methods are:

\begin{itemize}
\item {\bf "PCA"}: dimensionality reduction via principal component
analysis; the number of prinicipal components is determined such that at
least 95% of total variance in feature space can be explained 
\item {\bf "clustering"}: EM-clustering in feature space
\end{itemize}
}
\value{
If the function is called to automatically select prototypes, a character vector of gene IDs is returned.

If the function is called to perform dimensionality via PCA, the result is a list with items
\item{"features"}{ original data projected on the first k principal components }
\item{"pcs"}{ l x k matrix of principal components. Each column is one principal component }
\item{"lambda"}{ first k eigenvalues }

If the function is called to perform clustering in feature space, the
cluster centers are returned in a l x k matrix (each column is one
cluster center). The "Mclust" function in the package "mclust" is called
to perform the clustering. The BIC is used to calculate the optimal
number of clusters in the range 2,...,n. 
}
\references{
[1] H. Froehlich, N. Speer, C. Spieth, and A. Zell, Kernel Based Functional Gene Grouping, Proc. Int. Joint Conf. on Neural Networks (IJCNN), pp. 6886 - 6891, 2006\\
[2] N. Speer, H. Froehlich, A. Zell, Functional Grouping of Genes Using Spectral Clustering and Gene Ontology, Proc. Int. Joint Conf. on Neural Networks (IJCNN), pp. 298 - 303, 2005 }
\author{ Holger Froehlich }
\note{ The result depends on the currently set ontology ("BP","MF","CC").}
\seealso{ \code{\link{getGeneFeaturesPrototypes}}, \code{\link{getGeneSimPrototypes}}, \code{\link{setOntology}}  }
\examples{
\dontrun{ # takes too much time in the R CMD check
 proto=selectPrototypes(n=50) # --> returns a character vector of 50 genes with the highest number of annotations 
 feat=getGeneFeaturesPrototypes(c("207","208","7494"),prototypes=proto,pca=FALSE) # --> compute feature vectors 
 selectPrototypes(data=feat$features,method="pca") # ... and PCA projection
 }
}
\keyword{ file }% at least one, from doc/KEYWORDS
