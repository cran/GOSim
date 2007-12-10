\name{clusterEvaluation}
\alias{evaluateClustering}
\alias{analyzeCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evaluate a given grouping of genes or GO terms.  }
\description{
  \item{1.}{ Evaluate a given grouping of genes or terms with respect to their GO similarity.}
  \item{2.}{ Perform a GO enrichment analysis of a specific cluster using topGO }
}
\usage{
evaluateClustering(clust, Sim)

analyzeCluster(genesInCluster, allgenes, cutoff=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust}{ grouping = vector of integer or character }
  \item{Sim}{ similarity matrix }

  \item{genesInCluster}{Entrez gene IDs of genes in the cluster to investigate}
  \item{allgenes}{Entrez gene IDs of all genes (on the chip)}
  \item{cutoff}{significance cutoff for GO enrichment analysis}
}
\details{
   If necessary, more details than the description above 
}
\value{
  evaluateClustering returns a list with two items:
  \item{clusterstats }{matrix (ncluster x 2) of median within cluster similarities and median absolute deviations}
  \item{clustersil}{cluster silhouette values}

  analyzeCluster returns a list with two items:  
  \item{GOTerms}{list of significant GO terms and their description}
  \item{p.values}{vector of p-values for significant GO terms}
  \item{genes}{list of genes associated to each GO term}
}
\references{Rousseeuw, P., Silhouettes: a graphical aid to the interpretation and validation of cluster analysis, J. Comp. and Applied Mathematics, 1987, 20, 53-6

Adrian Alexa, J\"org Rahnenf\"uhrer, Thomas Lengauer: Improved scoring of functional groups from gene expression data by decorrelating GO graph structure, Bioinformatics, 2006, 22(13):1600-1607}
\author{Holger Froehlich}
\seealso{ \code{\link{getGeneSimPrototypes}}, \code{\link{getGeneSim}}, \code{\link{getTermSim}} }
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

	ev = evaluateClustering(cl, sim) # evaluate the clustering
	print(ev$clusterstats) # print out some statistics
	plot(ev$clustersil,main="") # plot the cluster silhouettes
	
	# investigate cluster 1 further	
	analyzeCluster(genesOfInterest[cl == 1], allgenes, cutoff=0.05) # print out what cluster 1 is about
	}
}
\keyword{ file }% at least one, from doc/KEYWORDS
