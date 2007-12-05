# evaluate a given clustering of genes or terms (e.g. terms significantly overrepresented in certain clusters of genes)  by means of the GO gene or term similarities
evaluateClustering<-function(clust, Sim){
	require(cluster)
	clus<-unique(clust)
	ncl<-length(clus)
	cluststats<-matrix(0,nrow=ncl,ncol=2)
	rownames(cluststats)<-clus
	colnames(cluststats)<-c("median within cluster similarity", "similarity mad")	
	for(c in 1:ncl){
		cl<-which(clust == clus[c])		
		S<-Sim[cl,cl]								
		cluststats[c,1]<-median(as.vector(S))
		cluststats[c,2]<-mad(as.vector(S))
	}
	sil<-silhouette(clust, as.dist(1-Sim))	
	list(clusterstats=cluststats,clustersil=sil)
}

# perform a GO gene set enrichment analysis for a specific cluster using topGO
analyzeCluster = function(genesInCluster, allgenes, cutoff=0.01){
	require(topGO)
	ontology = get("ontology", envir=GOSimEnv)	
	gomap <- get("gomap",env=GOSimEnv)
	anno <- gomap[as.character(allgenes)]
	goterms<-sapply(anno,function(x) names(x))
	geneList <- factor(as.integer(allgenes %in% genesInCluster))
	names(geneList) <- allgenes
	GOdata = new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = goterms)  
	test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.01) 
	res = getSigGroups(GOdata, test.stat) 
	sigterms = score(res)[score(res) < cutoff]
	goids<-as.list(GOTERM)
	require(annotate)
	goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]
	list(GOTerms=goids[names(sigterms)], p.values = sigterms)
}