#  perform PCA on data matrix X
# n<-number of PCs to extract (if not provided, use as many PCs as explain 95% of overall variance)
pca<-function(X, n=NULL){
  res<-eigen(cov(X), symmetric<-TRUE)
  lam<-res$values
  V<-res$vectors
  if(is.null(n)){
    n<-min(which(cumsum(lam)/sum(lam) >= 0.95))
    res<-sort(lam, index.return=TRUE, decreasing=TRUE)
  }
  V<-V[,res$ix[1:n]]
  return(list(features=X%*%V,pcs=V,lambda=lam[1:n]))
}
