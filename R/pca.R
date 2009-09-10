#  perform PCA on data matrix X
# n<-number of PCs to extract (if not provided, use as many PCs as explain 95% of overall variance)
pca<-function(X, n=NULL){	
	res = prcomp(X, center=TRUE, retx=TRUE)	
	lam = res$sdev^2
	if(is.null(n)){		
		n <- min(which(cumsum(lam)/sum(lam) >= 0.95))	
		n = max(n, 2)	
	}	
	s<-sort(lam, index.return=TRUE, decreasing=TRUE)
	features = res$x[,s$ix[1:n]]
	lambda = lam[s$ix[1:n]]
		
	list(features=features, pcs=res$rotation, lambda=lam)
}
