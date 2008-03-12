\name{setEnrichmentFactors}
\alias{setEnrichmentFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Set the depth and densitiy enrichment factors for GO term
  similarity. }
\description{
   Sets the depth and density enrichment factors for the enriched
   FuSSiMeg GO term similarity measure by Couto et al. 
}
\usage{
setEnrichmentFactors(alpha = 0.5, beta = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{ depth factor }
  \item{beta}{ density factor }
}
\value{
	none
}
\references{ F.Couto,M. Silva, P. Coutinho, Implementation of a
  Functional Semantic Similarity Measure between Gene-Products, DI/FCUL
  TR 03-29, Department of Informatics, University of Lisbon, 2003 } 
\author{ Holger Froehlich }
\note{ The enrichment factors are stored internally and are used by the
  function \code{\link{getTermSim}}, if one uses the method
  "CoutoEnriched" to calculate GO term similarities  
}
\seealso{ \code{\link{getTermSim}} }
\examples{
 if(require(graph)){
 	setEnrichmentFactors(alpha=0.1,beta=0.5)
 	getTermSim(c("GO:0006955","GO:0007584"),method="CoutoEnriched")
 }
}
\keyword{ file }% at least one, from doc/KEYWORDS

