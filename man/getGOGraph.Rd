\name{getGOGraph}
\alias{getGOGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get graphNEL object with leaves specified in the arguments.}
\description{
Returns a graphNEL object representing the GO graph with leaves specified in the argument. 
}
\usage{
getGOGraph(term)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{term}{ character vector of GO terms }
}
\details{
   The result is computed within the currently set ontology ("BP","MF","CC").  
}
\value{
  graphNEL object
}
\author{ Holger Froehlich }
\note{ directly calls the function \code{GOGraph} in the "GOstats" library 
}
\examples{
 G=getGOGraph(c("GO:0006955","GO:0007584"))
 library(Rgraphviz)
 plot(G)
}
\keyword{ file }% at least one, from doc/KEYWORDS
