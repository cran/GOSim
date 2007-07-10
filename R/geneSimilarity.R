# get GO information for a list of genes
getGOInfo<-function(geneIDs){	
	require(GO)
	if(!exists("GOSimEnv")) initialize()	
	ontology<-get("ontology",env=GOSimEnv)
	gomap<-get("gomap",env=GOSimEnv)
	goids<-as.list(GOTERM)
	require(annotate)
	goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]	
	anno<-gomap[as.character(geneIDs)]
	goterms<-sapply(anno,function(x) names(x))		
	if(class(goterms)=="list"){		
		goterms<- sapply(goterms,function(x) intersect(names(goids),x))
		info<-sapply(goterms, function(x) goids[x])
	}
	else{		
		goterms<-intersect(goterms, names(goids))
		info<-goids[goterms]
	}	
	info
}

# filter out genes not mapping to the category in question	
filterGO<-function(genelist){	
	if(!exists("GOSimEnv")) initialize()	
	IC<-get("IC", envir=GOSimEnv)
	ids<-names(IC[IC != Inf]) # only consider GO terms with some known annotation   
	gomap<-get("gomap",env=GOSimEnv)
	k<-1
	allgenes<-list()
	for(i in 1:length(genelist)){		
		annoi<-gomap[[as.character(genelist[i])]]						
		annoi<-intersect(names(annoi), as.character(ids))				
		if(length(annoi) > 0){
			allgenes[[k]]<-list(annotation=annoi,genename=as.character(genelist[i]))
			k<-k + 1
		}
	}	
	allgenes
}

# precompute term similarities for all pairs of GO terms belonging to the annotated genelists x (and y)
precomputeTermSims<-function(x, y=NULL, similarityTerm="JiangConrath", verbose=TRUE){ 	
	if(verbose)
		print("precomputing term similarities ...")	
	gotermsx<-unique(unlist(sapply(x, function(xx) xx$annotation)))
	if(!is.null(y)){      
		gotermsy<-unique(unlist(sapply(y, function(xx) xx$annotation)))		
		STerm<-matrix(0, nrow=length(gotermsx), ncol=length(gotermsy))
		rownames(STerm)=gotermsx
		colnames(STerm)=gotermsy
		for(i in 1:length(gotermsx)){						
			for(j in 1:length(gotermsy)){          
			STerm[i,j]<-calcTermSim(gotermsx[i],gotermsy[j], similarityTerm, verbose)
			}			
		}		
	} else{
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

# compute gene similarity for a pair of genes having GO terms anno1 and anno2
getGSim<-function(anno1, anno2, similarity="max", similarityTerm="JiangConrath", STerm=NULL, verbose=TRUE){	  
  if(length(anno1) < length(anno2)){
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
	# calculate term similarity
	ker<-matrix(0,nrow=length(a1),ncol=length(a2))	
	for(i in 1:length(a1)){
		for(j in 1:length(a2))
			ker[i,j]<-calcTermSim(a1[i],a2[j], similarityTerm, verbose)		
	}
  }  
  if(length(a1)*length(a2) > 0){
	if(similarity == "OA"){				
		res<-.C("OAWrapper", ker, nrow(ker), ncol(ker), as.integer(1), ret=double(1))$ret
		return(res)
	}
	else if(similarity == "max"){				
		return(max(ker))
	}
	else if(similarity == "mean"){				
		return(mean(ker))
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
getGeneSim<-function(genelist, similarity="OA", similarityTerm="JiangConrath", normalization=TRUE, verbose=TRUE){
	genelist <- unique(genelist)
	if(length(genelist) < 2)
		stop("Gene list should contain more than 2 elements!")
	allgenes<-filterGO(genelist)	
	if(length(allgenes) > 1){	
		STerm<-precomputeTermSims(x=allgenes, similarityTerm=similarityTerm, verbose=verbose) # precompute term similarities => speed up!		
		if(verbose)
			print(paste("Calculating similarity matrix with similarity measure",similarity))
		Ker<-matrix(0,nrow=length(allgenes),ncol=length(allgenes))
		colnames(Ker)<-sapply(allgenes,function(x) x$genename)
		rownames(Ker)<-colnames(Ker)
		for(i in 1:length(allgenes)){
			annoi<-(allgenes[[i]]$annotation)		
			Ker[i,i]<-getGSim(annoi,annoi, similarity, similarityTerm, STerm=STerm, verbose)
			if(i > 1){
				for(j in 1:(i-1)){
					annoj<-(allgenes[[j]]$annotation)
	# 			        print(paste(allgenes[[i]]$genename,allgenes[[j]]$genename))
					Ker[i,j]<-getGSim(annoi,annoj, similarity, similarityTerm, STerm=STerm, verbose)
					Ker[j,i]<-Ker[i,j]
				}
			}
		}			
		if(normalization){
			Kd<-sqrt(diag(Ker))
			Ker<-Ker/(Kd%*%t(Kd))
		}			
	}
	else{
		if(length(allgenes) == 0)
			stop("No gene has GO information!")					
		else if(length(allgenes) == 1)
			stop(paste("Only gene",allgenes," has GO information!"))					
	}
	Ker
}

# compute GO gene feature representation
getGeneFeaturesPrototypes<-function(genelist, prototypes=NULL, similarity="max", similarityTerm="JiangConrath", pca=TRUE, normalization=TRUE, verbose=TRUE){
	genelist<-unique(genelist)
	if(is.null(prototypes))
		prototypes<-selectPrototypes(verbose=verbose)
	allgenes<-filterGO(genelist)	
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
	if(pca){    
		pcares<-selectPrototypes(method="pca",data=PHI,verbose=verbose)		
		PHI<-pcares$features
		proto<-pcares$pcs    
	}
	if(normalization)
		PHI<-t(apply(PHI,1,function(x) return(x/(norm(x)+1e-10))))		  
	list(features=PHI,prototypes=proto)
}

# compute GO gene similarity using the feature representation
getGeneSimPrototypes<-function(genelist, prototypes=NULL, similarity="max", similarityTerm="JiangConrath", pca=TRUE, normalization=TRUE, verbose=TRUE){
	genelist<-unique(genelist)	
	if(length(genelist) < 2)
		stop("Gene list should contain more than 2 elements!")
	res<-getGeneFeaturesPrototypes(genelist, prototypes, similarity, similarityTerm, pca, normalization, verbose)
	PHI<-res$features		
	Ker<-PHI%*%t(PHI)
	if(normalization)
		Ker<-(Ker + 1)/2 # scale to [0,1]
	list(similarity=Ker,prototypes=res$prototypes,features=res$features)
}

# method to a) select prototype genes b) to perform subselection of prototypes via i) PCA ii) clustering
selectPrototypes<-function(n=250, method="frequency", data=NULL, verbose=TRUE){
	if(!exists("GOSimEnv")) initialize()	
	ontology<-get("ontology",env=GOSimEnv)  
	if(method == "frequency"){
		if(verbose)
			print("Automatic determination of prototypes using most frequent genes ...")
		locus<-as.list(GOENTREZID)
		goids<-as.list(GOTERM)
		goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]		
		locusnames<-intersect(names(locus),names(goids))
		locus<-locus[locusnames]
		lt<-table(unlist(locus))
		res<-sort(lt,decreasing=TRUE)
		prototypes<-names(res)[1:n]
		return(prototypes)
	}
	else if(method == "random"){
		if(verbose)		
			print("Automatic determination of prototypes using random genes ...")
		locus<-as.list(GOENTREZID)
		goids<-as.list(GOTERM)
		goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]		
		locusnames<-intersect(names(locus),names(goids))
		locus<-locus[locusnames]
		lt<-names(table(unlist(locus)))
		set.seed(0)
		ridx<-round(runif(n,min=1,max=length(lt)))
		prototypes<-lt[ridx]		
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
		require(mclust)
		if(is.null(data))
			stop("You need to specify a data matrix with feature vectors")
		if(verbose)
			print("Clustering feature vectors ...")
		res<-Mclust(t(data),2,n)
		return(res$mu)
	}
	else
	stop(paste("selectPrototypes: Unknown method",method,"!"))
}

