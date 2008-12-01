\name{clusterEvaluation}
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
  \item{clust}{ vector of cluster labels (integer or character) for each gene }
  \item{Sim}{ similarity matrix }
}
\details{
   If necessary, more details than the description above 
}
\value{
  evaluateClustering returns a list with two items:
  \item{clusterstats }{matrix (ncluster x 2) of median within cluster similarities and median absolute deviations}
  \item{clustersil}{cluster silhouette values}  
}
\references{Rousseeuw, P., Silhouettes: a graphical aid to the interpretation and validation of cluster analysis, J. Comp. and Applied Mathematics, 1987, 20, 53-6
}
\author{Holger Froehlich}
\seealso{ \code{\link{getGeneSimPrototypes}}, \code{\link{getGeneSim}}, \code{\link{getTermSim}}, \code{\link{GOenrichment}} }
\examples{
	\dontrun{	
	setOntology("BP")
	gomap <- get("gomap",env=GOSimEnv)
	allgenes = sample(names(gomap), 1000) # suppose these are all genes
	genesOfInterest = sample(allgenes, 20) # suppose these are all genes of interest
	
	sim = getGeneSim(genesOfInterest,verbose=FALSE) # and these are their similarities
	hc = hclust(as.dist(1-sim), method="ward") # use them to perform a clustering
	plot(hc) # plot the cluster tree
	cl = cutree(hc, k=3) # take 3 clusters	

	if(require(cluster)){
		ev = evaluateClustering(cl, sim) # evaluate the clustering
		print(ev$clusterstats) # print out some statistics
		plot(ev$clustersil,main="") # plot the cluster silhouettes
	}
	
	# investigate cluster 1 further	
	if(require(topGO))
		GOenrichment(genesOfInterest[cl == 1], allgenes, cutoff=0.05) # print out what cluster 1 is about
	}
}
\keyword{ file }% at least one, from doc/KEYWORDS
