\name{getParents}
\alias{getParents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get direct parents for each GO term. }
\description{
   Returns the list of all direct parents ( = less specific terms one hiearchy level up) associated to each GO term. The type of relationship between GO terms ("is a" or "part of") is ignored.
}
\usage{
getParents()
}
%- maybe also 'usage' for other objects documented here.
\value{
List with entry names for each GO term. Each entry contains a character vector with the direct parent GO terms. 
}
\author{ Holger Froehlich }
\note{ The result is computed within the currently set ontology
  ("BP","MF","CC"). It directly uses the "GO" library to compute the
  result. 
}
\seealso{ \code{\link{getOffsprings}}, \code{\link{getChildren}}, \code{\link{getAncestors}}, \code{\link{setOntology}} }
\examples{
 \dontrun{getParents()}
}
\keyword{ file }% at least one, from doc/KEYWORDS
