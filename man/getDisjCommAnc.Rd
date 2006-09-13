\name{getDisjCommAnc}
\alias{getDisjCommAnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get disjoint common ancestors. }
\description{
   Returns the GO terms representing the disjoint common ancestors of two GO terms. 
}
\usage{
getDisjCommAnc(term1, term2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{term1}{ GO term 1 }
  \item{term2}{ GO term 2 }
}
\details{
   The result is computed within the currently set ontology ("BP","MF","CC"). 
}
\value{
Character vector of GO terms
}
\references{ Couto, F.; Silva, M. & Coutinho, P., Semantic Similarity over the Gene Ontology: Family Correlation and Selecting Disjunctive Ancestors, Conference in Information and Knowledge Management, 2005 }
\author{ Holger Froehlich }

\seealso{ \code{\link{getTermSim}},\code{\link{getGOGraph}}, \code{link{setOntology}}  }
\examples{
	getDisjCommAnc("GO:0006955","GO:0007584")
}
\keyword{ file }% at least one, from doc/KEYWORDS
