\name{getOffsprings}
\alias{getOffsprings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get all offspring associated with one or more GO term }
\description{
   Returns the list of all (also indirect) offspring (= more specific terms) associated to each GO term. The type of relationship between GO terms ("is a" or "part of") is ignored.
}
\usage{
getOffsprings()
}
%- maybe also 'usage' for other objects documented here.
\value{
List with entry names for each GO term. Each entry contains a character
vector with the offspring GO terms.
}
\author{ Holger Froehlich }
\note{ The result is computed within the currently set ontology
  ("BP","MF","CC"). It directly uses the "GO" library to compute the
  result. 
}
\seealso{ \code{\link{getChildren}}, \code{\link{getParents}}, \code{\link{getAncestors}}, \code{\link{setOntology}} }
\examples{
 \dontrun{getOffsprings()}
}
\keyword{ file }% at least one, from doc/KEYWORDS
