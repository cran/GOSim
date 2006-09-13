\name{evaluateClustering}
\alias{evaluateClustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evaluate a given grouping of genes or GO terms.  }
\description{
  Evaluate a given grouping of genes or terms with respect to their GO similarity.
}
\usage{
evaluateClustering(clust, Sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust}{ grouping = vector of integer or character }
  \item{Sim}{ similarity matrix }
}
\details{
   If necessary, more details than the description above 
}
\value{
  List with two items:
  \item{clusterstats }{matrix (ncluster x 2) of median within cluster similarities and median absolute deviations}
  \item{clustersil}{cluster silhouette values}  
}
\references{Rousseeuw, P., Silhouettes: a graphical aid to the interpretation and validation of cluster analysis, J. Comp. and Applied Mathematics, 1987, 20, 53-6}
\author{Holger Froehlich}
\seealso{ \code{\link{getGeneSimPrototypes}}, \code{\link{getGeneSim}}, \code{\link{getTermSim}} }
\examples{
	setOntology("BP")
	genes = c("207","208","596","901","780","3169","9518","2852","26353","8614","7494")
	sim = getGeneSim(genes,verbose=FALSE)
	ev = evaluateClustering(c(2,3,2,3,1,2,1,1,3,1,2), sim)
	print(ev$clusterstats)
	plot(ev$clustersil,main="")
}
\keyword{ file }% at least one, from doc/KEYWORDS
