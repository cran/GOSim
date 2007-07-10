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
