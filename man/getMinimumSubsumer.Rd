\name{getMinimumSubsumer}
\alias{getMinimumSubsumer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute minimum subsumer of two GO terms.}
\description{
Returns the minimum subsumer (i.e. the common ancestor having the
maximal information content) of two GO terms 
}
\usage{
getMinimumSubsumer(term1, term2)
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
  GO term representing the minimum subsumer
}
\references{P. Resnik, Using Information Content to evaluate semantic similarity in a taxonomy, Proc. 14th Int. Conf. Artificial Intel., 1995 }
\author{ Holger Froehlich}
\seealso{ \code{\link{getTermSim}}, \code{\link{getGOGraph}}, \code{link{setOntology}} }
\examples{
 setOntology("BP")
 getMinimumSubsumer("GO:0006955","GO:0007584")
 # returns GO:0050896
}
\keyword{ file }% at least one, from doc/KEYWORDS
