\name{getChildren}
\alias{getChildren}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get a list of all direct children of each GO term. }
\description{
   Returns the list of all direct children (= more specific terms one hierarchy level down)
   associated to each GO term. The type of relationship between GO terms ("is a" or "part of") is ignored.
}
\usage{
getChildren()
}
%- maybe also 'usage' for other objects documented here.
\value{
List with entry names for each GO term. Each entry contains a character
vector with the direct children GO terms. 
}
\author{ Holger Froehlich }
\note{ The result is computed within the currently set ontology
  ("BP","MF","CC"). It directly uses the "GO" library to compute the
  result.
}
\seealso{ \code{\link{getOffsprings}}, \code{\link{getParents}},
  \code{\link{getAncestors}}, \code{\link{setOntology}} }
\examples{
	\dontrun{getChildren()}
}
\keyword{ file }% at least one, from doc/KEYWORDS
