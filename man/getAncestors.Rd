\name{getAncestors}
\alias{getAncestors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ get list of ALL ancestors associated to each GO term }
\description{
   Returns the list of all (also indirect) ancestors (= less specific terms) associated to each GO term. The type of relationship between GO terms ("is a" or "part of") is ignored.
}
\usage{
getAncestors()
}
%- maybe also 'usage' for other objects documented here.
\value{
List with entry names for each GO term. Each entry contains a character
vector with the ancestor GO terms. 
}
\author{ Holger Froehlich }
\note{ The result is computed within the currently set ontology
  ("BP","MF","CC"). It directly uses the "GO" library to compute the
  result. 
}
\seealso{ \code{\link{getOffsprings}}, \code{\link{getChildren}}, \code{\link{getParents}}, \code{\link{setOntology}} }
\examples{
	\dontrun{getAncestors()}
}
\keyword{ file }% at least one, from doc/KEYWORDS
