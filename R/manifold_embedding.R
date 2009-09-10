calc.diffusion.kernel = function(method="diffKernelLapl", m=7){		
	require(GO.db)
	require(Matrix)	
	require(corpcor)
	if(!exists("GOSimEnv")) initialize()
	evidences<-get("evidences", envir=GOSimEnv)
	ontology<-get("ontology",envir=GOSimEnv)
	organism = get("organism", envir=GOSimEnv)
	print(paste("calculating diffusion kernel '", method, "' for ontology ", ontology, " using evidence codes '", paste(evidences,collapse=", "), "' (", organism,") ...",sep=""))		
	ids<- toTable(GOTERM)	 
	ids = unique(ids[ids[,"Ontology"] == ontology,"go_id"]) # these are all GO terms, which belong to the corrseponding category	
	# GO term frequency: freq(subgraph rooted at t): each gene function is viewed as ONE GO graph
	ancestor = getAncestors()	
	gomap <- get("gomap",env=GOSimEnv)			
	goterms = sapply(gomap, function(g){
		myterms = names(g)
		myterms = intersect(myterms, ids)
		unique(c(myterms, unlist(ancestor[myterms])))
	})
	tab = table(unlist(goterms))	
	# GO term frequency: freq(t)
	goterms2 <- unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) # all GO terms appearing in an annotation
	goterms2 <- goterms2[goterms2 %in% ids] # this is to ensure that we only get GO terms mapping to the given ontology
	tab2 <- table(goterms2)	
	G = getGOGraph(intersect(names(tab),names(tab2)), Inf) # edges go from specific to more general terms!
	W = as(G, "matrix")	
	W = W[1:(nrow(W)-1),1:(ncol(W)-1)] # cut row/column "all"	
	edges = which(W == 1, arr.ind=TRUE)	
	for(i in 1:nrow(edges)){
		if(ids[edges[i,2]] %in% names(tab2))
			freq_a = tab2[ids[edges[i,2]]]
		else
			freq_a = 0
		if(ids[edges[i,1]] %in% names(tab))
			freq_B = tab[ids[edges[i,1]]]
		else
			freq_B = 0
		W[edges[i,1],edges[i,2]] = freq_a + freq_B + 1# freq(a->b) = freq(a) + freq(subgraph rooted at b)
		W[edges[i,2],edges[i,1]] = W[edges[i,1],edges[i,2]]
	}	
	W[W != 0] = 1/W[W != 0]
	W = W/max(W)
	W = Matrix(W, sparse=TRUE)
	deg = rowSums(W)
	D = Diagonal(x=deg)
	if(method == "diffKernelgraphLapl")
		K = pseudoinverse(D - W)
	else if(method == "diffKernelpower"){
		invD_half = Diagonal(x=1/sqrt(deg))
		P2 = invD_half %*% W %*% invD_half
		K = P2%*%P2
		if(m > 1){
			for(i in 2:m)
				K = K%*%P2%*%P2
		}		
	}
	else if(method %in% c("diffKernelLLE")){
		invD = Diagonal(x=1/deg)
		P = invD%*%W
		I = Diagonal(ncol(P)) 
		T = I - P
		M = t(T)%*%(T)	
		lam = eigen(M, symmetric=TRUE, only.values=TRUE)
		K = lam$values[1]*I - M
		E = I - matrix(1,ncol=ncol(I), nrow=nrow(I))/ncol(I)
		K = E%*%K%*%E						
	}
	else if(method == "diffKernelexpm"){
		L = D - W
		K = expm(-m*L)
	}
	K = 0.5*(K + t(K)) # force it to be symmetric
	K = normalize.kernel(K, "sqrt")
	K = (K + 1)/2		
	dimnames(K) = dimnames(W)
	save(K, file=paste(method, ontology,organism,paste(evidences,collapse="_"),".rda",sep=""))
	print("done.")	
}

load.diffusion.kernel = function(method="diffKernelLapl"){
	evidences<-get("evidences", envir=GOSimEnv)
	ontology<-get("ontology",envir=GOSimEnv)
	organism = get("organism", envir=GOSimEnv)
	fname = paste(method,ontology,organism, paste(evidences,collapse="_"),sep="")	
	tryCatch(utils::data(list=fname,package="GOSim",envir=GOSimEnv), warning=function(w) stop(paste("File", fname, "with diffusion kernel", method, "for organism '", organism, "', ontology '", ontology, "', evidence codes '", evidences, "' not found!\nPlease invoke calc.diffusion.kernel() to calculate the diffusion kernel!")))
}
