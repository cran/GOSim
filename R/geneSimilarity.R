# get GO information for a list of genes
getGOInfo<-function(geneIDs){	
	geneIDs = as.character(geneIDs)
	require(GO.db)
	if(!require(annotate))
		stop("Package annotate is required for function getGOInfo")
	if(!exists("GOSimEnv")) initialize()		
	ontology<-get("ontology",env=GOSimEnv)
	gomap<-get("gomap",env=GOSimEnv)	
	gomap = gomap[geneIDs]		
	goterms = lapply(gomap, function(x) sapply(x, function(y) y$Ontology == ontology))
	goterms = goterms[!is.na(names(goterms))]
	if(length(goterms) == 0)
		stop("No GO information available for these genes!")
	goterms <-sapply(goterms, function(x) unique(names(x[which(x)])))
	goterms <- goterms[sapply(goterms,length) > 0]
	goids<- toTable(GOTERM)				
	IC<-get("IC", envir=GOSimEnv)	
	info = sapply(goterms, function(g) cbind(unique(goids[goids[,"go_id"] %in% g,c("go_id", "Term", "Definition")]), IC=IC[g]))
	info
}

# filter out genes not mapping to the category in question	
filterGO<-function(genelist){	
	cat("filtering out genes not mapping to the currently set GO category ...")
	if(!exists("GOSimEnv")) initialize()	
	IC<-get("IC", envir=GOSimEnv)
	ids<-names(IC[IC != Inf]) # only consider GO terms with some known annotation   
	gomap<-get("gomap",env=GOSimEnv)
	k<-1
	allgenes<-list()
	for(i in 1:length(genelist)){		
		annoi<-gomap[[match(as.character(genelist[i]), names(gomap))]]
		annoi<-intersect(names(annoi), as.character(ids))				
		if(length(annoi) > 0){
			allgenes[[k]]<-list(annotation=annoi,genename=as.character(genelist[i]))
			k<-k + 1
		}
	}	
	cat(" ===> list of ", length(genelist), "genes reduced to ", length(allgenes), "\n")
	allgenes
}

#calc.sim.significance = function(Sim, B=1000, similarity="funSimMax", similarityTerm="relevance", normalization=TRUE, method="sqrt", avg=(similarity=="OA"), adj.method="bonf"){	
#	n = NCOL(Sim)
#	pvals = matrix(1, ncol=n, nrow=n)			
#	gomap <- get("gomap",env=GOSimEnv)
#	allgenes = filterGO(names(gomap))
#	allgenes = sapply(allgenes, function(x) x$genename)		
#	for(b in 1:B){
#		sam = base:::sample(allgenes, n)
#		K = getGeneSim(sam, similarity=similarity, similarityTerm=similarityTerm, normalization=normalization, method=method, avg=avg, verbose=FALSE)
#		pvals = (K >= Sim)*1 + pvals					
#	}		
#	pvals = pvals/B
#	pvals = p.adjust(pvals, method=adj.method)
#	pvals
#}

# get GO graphs for GO terms associated to each gene
getGOGraphsGenes <- function(genelist, prune=Inf){
	allgenes = filterGO(genelist)
	if(length(allgenes) > 0){
		G = list()
		for(i in 1:length(allgenes)){
			G[[i]] = getGOGraph(allgenes[[i]]$annotation, prune)
		}
	}
	else
		stop("No gene has GO information!")
	G
}

# precompute term similarities for all pairs of GO terms belonging to the annotated genelists x (and y)
precomputeTermSims<-function(x, y=NULL, similarityTerm="JiangConrath", verbose=FALSE){ 	
	if(verbose)
		print("precomputing term similarities ...")	
	gotermsx<-as.vector(unique(unlist(sapply(x, function(xx) xx$annotation))))
	if(!is.null(y)){  		 
		gotermsy<-as.vector(unique(unlist(sapply(y, function(xx) xx$annotation))))
		if(similarityTerm %in% c("diffKernelgraphLapl", "diffKernelLLE", "diffKernelpower", "diffKernelexpm")){
			STerm = getTermSim(c(gotermsx, gotermsy), method=similarityTerm)
			STerm = STerm[gotermsx, gotermsy]
			return(STerm)
		}
		STerm<-matrix(0, nrow=length(gotermsx), ncol=length(gotermsy))
		rownames(STerm)=gotermsx
		colnames(STerm)=gotermsy
		for(i in 1:length(gotermsx)){						
			for(j in 1:length(gotermsy)){          
				STerm[i,j]<-calcTermSim(gotermsx[i],gotermsy[j], similarityTerm, verbose)
			}			
		}		
	} else{
		if(similarityTerm %in% c("diffKernelgraphLapl", "diffKernelLLE", "diffKernelpower", "diffKernelexpm")){
			STerm = getTermSim(gotermsx, method=similarityTerm)			
			return(STerm)
		}
		STerm<-matrix(0, nrow=length(gotermsx), ncol=length(gotermsx))
		rownames(STerm)<-gotermsx
		colnames(STerm)<-gotermsx
		for(i in 1:length(gotermsx)){
			STerm[i,i]<-calcTermSim(gotermsx[i], gotermsx[i], similarityTerm, verbose)
			if(i > 1){
				for(j in 1:(i-1)){					
					STerm[i,j]<-calcTermSim(gotermsx[i], gotermsx[j], similarityTerm, verbose)
					STerm[j,i]<-STerm[i,j]
				}
			}
		}
	}	
	STerm
}

getWeightedDotSim <- function(anno1, anno2){
	v1 = getGeneFeatures.internal(anno1)
	v2 = getGeneFeatures.internal(anno2)
	dot = crossprod(v1,v2)	
}

# compute gene similarity for a pair of genes having GO terms anno1 and anno2
getGSim<-function(anno1, anno2, similarity="max", similarityTerm="JiangConrath", STerm=NULL, avg=FALSE, verbose=FALSE){	  
  if(length(anno1) <= length(anno2)){
	a1<-anno1
	a2<-anno2
	swap<-FALSE
  }
  else{
	a1<-anno2
	a2<-anno1
	swap<-TRUE
  }		
  if(!is.null(STerm)){ # use precomputed similarity values
	if(!swap)
		ker<-STerm[a1,a2]
	else
		ker<-STerm[a2,a1]		
	if(length(a1) == 1)
		ker<-t(as.matrix(ker))		
	if(is.null(ker) || is.null(nrow(ker))){
		warning(paste("No GO information for",a1,a2,". Similarity set to NaN."))		
		return(NaN)
	}	
	if(nrow(ker) > ncol(ker))
		ker<-t(ker)
  }
  else{ 
	if(similarity %in% c("dot"))
		return(getWeightedDotSim(a1, a2))			
	else{
		# calculate term similarity
		ker<-matrix(0,nrow=length(a1),ncol=length(a2))	
		for(i in 1:length(a1)){
			for(j in 1:length(a2))
				ker[i,j]<-calcTermSim(a1[i],a2[j], similarityTerm, verbose)		
		}
	}
  }  
  if(length(a1)*length(a2) > 0){
	if(similarity == "OA"){				
		res<-.C("OAWrapper", ker, nrow(ker), ncol(ker), as.integer(1), ret=double(1))$ret
		if(avg)
			res = res/length(a2)	
		return(res)
	}
	else if(similarity == "max"){				
		return(max(ker))
	}
	else if(similarity == "mean"){				
		return(mean(ker))
	}  
	else if(similarity == "funSimAvg"){
		rowMax = mean(apply(ker,1,max))
		colMax = mean(apply(ker,2,max))
		return(0.5*(rowMax + colMax))
	}
	else if(similarity == "funSimMax"){
		rowMax = mean(apply(ker,1,max))
		colMax = mean(apply(ker,2,max))
		return(max(rowMax, colMax))
	}
	else if(similarity == "hausdorff"){
		rowMax = min(apply(ker,1,max))
		colMax = min(apply(ker,2,max))
		return(min(rowMax, colMax))	
	}	
	else
		stop(paste("getGSim: Unknown gene similarity",similarity,"!"))
  }
  else{	
	warning(paste("No GO information for",a1,a2,". Similarity set to NaN."))		
	return(NaN)
  }
}

# compute GO gene similarity for a list of genes
# getGeneSim<-function(genelist, similarity="funSimMax", similarityTerm="Lin", normalization=!(similarity %in% c("funSimAvg","funSimMax")), method="sqrt", avg=(similarity=="OA"), verbose=FALSE){
getGeneSim<-function(genelist1, genelist2=NULL, similarity="funSimMax", similarityTerm="relevance", normalization=TRUE, method="sqrt", avg=(similarity=="OA"), verbose=FALSE){	
	genelist1 <- unique(genelist1)
	genelist2 = unique(genelist2)
	if(length(genelist1) < 2 && is.null(genelist2))
		stop("Gene list should contain more than 2 elements!")
	allgenes<-filterGO(genelist1)	
	if(length(allgenes) > 1 && is.null(genelist2)){			
		if(!(similarity %in% c("dot")))
			STerm<-precomputeTermSims(x=allgenes, similarityTerm=similarityTerm, verbose=verbose) # precompute term similarities => speed up!		
		else
			STerm = NULL
		if(verbose)
			print(paste("Calculating similarity matrix with similarity measure",similarity))
		Ker<-matrix(0,nrow=length(allgenes),ncol=length(allgenes))
		colnames(Ker)<-sapply(allgenes,function(x) x$genename)
		rownames(Ker)<-colnames(Ker)
		for(i in 1:length(allgenes)){
			annoi<-(allgenes[[i]]$annotation)		
			Ker[i,i]<-getGSim(annoi,annoi, similarity, similarityTerm, STerm=STerm, avg=avg, verbose)
			if(i > 1){
				for(j in 1:(i-1)){
					annoj<-(allgenes[[j]]$annotation)
	# 			        print(paste(allgenes[[i]]$genename,allgenes[[j]]$genename))
					Ker[i,j]<-getGSim(annoi,annoj, similarity, similarityTerm, STerm=STerm, avg=avg, verbose)
					Ker[j,i]<-Ker[i,j]
				}
			}
		}			
		if(normalization){			
			Ker = normalize.kernel(Ker, method=method)
			if(any(Ker > 1, na.rm=T)) # this has been updated
				warning("Similarity matrix contains values > 1! This may happen with simlarity='funSimMax', if one gene's GO annotation is a complete subset of another gene's GO annotation.")
			Ker[Ker>1] = 1 # can happen with similarity funSimMax in cases where one GO annotation is subset of another one
		}			
	}
	else if(length(allgenes) > 0 && length(genelist2) > 0){
		allgenes2<-setdiff(filterGO(genelist2), allgenes)
		if(!(similarity %in% c("dot")))
			STerm<-precomputeTermSims(x=allgenes, y=allgenes2, similarityTerm=similarityTerm, verbose=verbose) # precompute term similarities => speed up!		
		else
			STerm = NULL
		if(verbose)
			print(paste("Calculating similarity matrix with similarity measure",similarity))
		Ker<-matrix(0,nrow=length(allgenes),ncol=length(allgenes2))
		colnames(Ker)<-sapply(allgenes2,function(x) x$genename)
		rownames(Ker)<-sapply(allgenes,function(x) x$genename)
		for(i in 1:length(allgenes)){
			annoi<-(allgenes[[i]]$annotation)
			if(normalization)
				kerselfi = getGSim(annoi, annoi, similarity, similarityTerm, STerm=NULL, avg=avg, verbose)
			for(j in 1:length(allgenes2)){
				annoj<-(allgenes2[[j]]$annotation)
				# 			        print(paste(allgenes[[i]]$genename,allgenes[[j]]$genename))
				Ker[i,j]<-getGSim(annoi,annoj, similarity, similarityTerm, STerm=STerm, avg=avg, verbose)	
				if(normalization){
					kerselfj = getGSim(annoj, annoj, similarity, similarityTerm, STerm=NULL, avg=avg, verbose)
					Ker[i,j] = normalize.kernel(Ker[i,j], kerselfi, kerselfj, method=method)
				}
			}
		}			
		if(any(Ker > 1, na.rm=T))
			warning("Similarity matrix contains values > 1! This may happen with simlarity='funSimMax', if one gene's GO annotation is a complete subset of another gene's GO annotation.")
		Ker[Ker>1] = 1 # can happen with similarity funSimMax in cases where one GO annotation is subset of another one
	}
	else{
		if(length(allgenes) == 0)
			stop("No gene has GO information!")					
		else if(length(allgenes) == 1)
			stop(paste("Only gene",allgenes," has GO information!"))					
	}
	Ker
}

getGeneFeatures.internal = function(anno){
	ancestor<-get("ancestor",envir=GOSimEnv)	
	an<-unlist(ancestor[names(ancestor) %in% anno])	
	IC<-get("IC", envir=GOSimEnv)
	v = double(length(IC))	
	names(v) = names(IC)	
	v[c(anno, an)] = IC[c(anno, an)]
	v
}

# compute feature representation for genes
getGeneFeatures = function(genelist, pca=FALSE, normalization=FALSE, verbose=FALSE){
	if(!exists("GOSimEnv")) initialize()	
	genelist <- unique(genelist)
	if(length(genelist) < 1)
		stop("Gene list should contain at least 1 element!")
	allgenes<-filterGO(genelist)	
	IC<-get("IC", envir=GOSimEnv)
	if(length(allgenes) > 0){			
		PHI<-matrix(0,nrow=length(allgenes),ncol=length(IC))
		rownames(PHI) <- sapply(allgenes,function(x) x$genename)
		colnames(PHI) <- names(IC)
		for(i in 1:length(allgenes)){
			annoi <- (allgenes[[i]]$annotation)		
			PHI[i,] <- getGeneFeatures.internal(annoi)							
		}	
		if(pca){
			pcares<-selectPrototypes(method="pca",data=PHI,verbose=verbose)		
			PHI<-pcares$features
		}				
		if(normalization)
			PHI<-t(apply(PHI,1,function(x) return(x/(norm(x)+1e-10))))					
	}
	else{		
		stop("No gene has GO information!")							
	}
	PHI
}

# compute prototype feature representation for genes
getGeneFeaturesPrototypes<-function(genelist, prototypes=NULL, similarity="max", similarityTerm="JiangConrath", pca=TRUE, normalization=TRUE, verbose=FALSE){
	genelist<-unique(genelist)
	if(is.null(prototypes))
		prototypes<-selectPrototypes(verbose=verbose)
	allgenes<-filterGO(genelist)
	if(length(allgenes) == 0)
		stop("No gene has GO information!")
	proto<-filterGO(prototypes)
	if(length(proto) == 0)
		stop("getGeneFeaturesPrototypes: Number of prototypes equals zero after filtering!")	
	STerm<-precomputeTermSims(x=allgenes,y=proto,similarityTerm=similarityTerm,verbose=verbose)  # precompute term similarities => speed up!
	if(verbose)
		print(paste("calculating feature vectors with",length(proto),"prototypes"))
	PHI<-matrix(0,nrow=length(allgenes),ncol=length(proto))
	for(i in 1:length(allgenes)){		
		for(j in 1:length(proto)){
			annoi<-(allgenes[[i]]$annotation)		
			annoj<-(proto[[j]]$annotation)			
			PHI[i,j]<-getGSim(annoi,annoj, similarity, similarityTerm, STerm=STerm, verbose)
		}
	}			
	rownames(PHI)<-sapply(allgenes,function(x) x$genename)
	colnames(PHI)<-sapply(proto,function(x) x$genename)	
	if(pca & (length(proto) > 2)){    
		pcares<-selectPrototypes(method="pca",data=PHI,verbose=verbose)		
		PHI<-pcares$features
		proto<-pcares$pcs    
	}	
	if(normalization)
		PHI<-t(apply(PHI,1,function(x) return(x/(norm(x)+1e-10))))		  
	list(features=PHI,prototypes=proto)
}

# compute GO gene similarity using the prototype feature representation
getGeneSimPrototypes<-function(genelist, prototypes=NULL, similarity="max", similarityTerm="JiangConrath", pca=TRUE, normalization=TRUE, verbose=FALSE){
	genelist<-unique(genelist)	
	if(length(genelist) < 2)
		stop("Gene list should contain at least 2 elements!")
	res<-getGeneFeaturesPrototypes(genelist, prototypes, similarity, similarityTerm, pca, normalization, verbose)
	PHI<-res$features		
	Ker<-PHI%*%t(PHI)
	if(normalization)
		Ker<-(Ker + 1)/2 # scale to [0,1]
	list(similarity=Ker,prototypes=res$prototypes,features=res$features)
}

# method to a) select prototype genes b) to perform subselection of prototypes via i) PCA ii) clustering
selectPrototypes<-function(n=250, method="frequency", data=NULL, verbose=FALSE){
	if(!exists("GOSimEnv")) initialize()	
	ontology<-get("ontology",env=GOSimEnv) 	
	if(method == "frequency"){
		if(verbose)
			print("Automatic determination of prototypes using genes with most frequent annotation ...")	
		gomap<-get("gomap",env=GOSimEnv)
		goterms = lapply(gomap, function(x) sapply(x, function(y) y$Ontology == ontology))
		goterms = goterms[!is.na(names(goterms))]
		if(length(goterms) == 0)
			stop("No GO information available for these genes!")
		goterms <-sapply(goterms, function(x) names(x[which(x)]))
		freq <- sapply(goterms,length)
		res = sort(freq, decreasing=TRUE)		
# 		locus<-as.list(GOENTREZID)
# 		goids<-as.list(GOTERM)
# 		goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]		
# 		locusnames<-intersect(names(locus),names(goids))
# 		locus<-locus[locusnames]
# 		lt<-table(unlist(locus))
# 		res<-sort(lt,decreasing=TRUE)
		prototypes<-names(res)[1:n]
		return(prototypes)
	}
	else if(method == "random"){
		if(verbose)		
			print("Automatic determination of prototypes using random genes from the current ontology...")
		gomap<-get("gomap",env=GOSimEnv)
		goterms = lapply(gomap, function(x) sapply(x, function(y) y$Ontology == ontology))
		goterms = goterms[!is.na(names(goterms))]
		if(length(goterms) == 0)
			stop("No GO information available for these genes!")
		goterms <-sapply(goterms, function(x) names(x[which(x)]))
		goterms <- goterms[sapply(goterms,length) > 0]
		prototypes = sample(names(goterms), n)
# 		locus<-as.list(GOENTREZID)
# 		goids<-as.list(GOTERM)
# 		goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]		
# 		locusnames<-intersect(names(locus),names(goids))
# 		locus<-locus[locusnames]
# 		lt<-names(table(unlist(locus)))
# 		set.seed(0)
# 		ridx<-round(runif(n,min=1,max=length(lt)))
# 		prototypes<-lt[ridx]		
		return(prototypes)
	}
	else if(method == "pca"){
		if(is.null(data))
			stop("You need to specify a data matrix with feature vectors")
		pcares<-pca(data)
		if(verbose)
			print(paste("PCA: dimension reduced to",ncol(pcares$features),"principal components (95% total variance explained)"))
		return(pcares)
	}
	else if(method == "clustering"){    		
		if(is.null(data))
			stop("You need to specify a data matrix with feature vectors")
		if(verbose)
			print("Clustering feature vectors ...")
		#res<-Mclust(t(data),2:n)
		#return(res$mu)
		res = sapply(2:n, function(k) flexmix(t(data)~1, cluster=cutree(hclust(dist(t(data))),k=k), k=k, model= FLXMCmvnorm()))
		bics = sapply(res, BIC)
		res = res[[which.min(bics)]]
		p = parameters(res)
		return(p[setdiff(1:NROW(p), grep("cov", rownames(p))),])
	}
	else
	stop(paste("selectPrototypes: Unknown method",method,"!"))
}

