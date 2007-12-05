\name{getGOGraph}
\alias{getGOGraph}
\alias{getGOGraphsGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get graphNEL object with leaves specified in the arguments.}
\description{
The function getGOGraph returns a graphNEL object representing the GO graph with leaves specified in the argument. The function getGOGraphsGenes returns a set of graphNEL objects. The ith graph object is created by call to getGOGraph with the GO terms associated to gene i. It hence shows for each gene, where its GO terms are located within the GO structure.
}
\usage{
getGOGraph(term, prune=Inf)

getGOGraphsGenes(genelist, prune=Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{term}{ character vector of GO terms }
  \item{genelist}{ character vector of Entrez gene IDs}
  \item{prune}{ do not show the complete graph, but prune it after the specified number of ancestors }
}
\details{
   The result is computed within the currently set ontology ("BP","MF","CC").  
}
\value{
  graphNEL object(s)
}
\author{ Holger Froehlich }
\note{ directly calls the function \code{GOGraph} in the "GOstats" library 
}
\examples{
\dontrun{
 G=getGOGraph(c("GO:0006955","GO:0007584"))
 library(Rgraphviz)
 plot(G)
 Gs = getGOGraphsGenes(c("207","7494"))
 plot(Gs[[1]]) # plot the first of both GO graphs
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
