setEnrichmentFactors<-function(alpha=0.5, beta=0.5){
	if(!exists("GOSimEnv")) initialize()
	assign("alphaParam",alpha, envir=GOSimEnv)
	assign("betaParam",beta, envir=GOSimEnv)
}

GOGraph = function(term, env){
    oldEdges <- vector("list", length = 0)
    oldNodes <- vector("character", length = 0)
    newN <- term
    done <- FALSE
    while (!done) {
        newN <- newN[!(newN %in% oldNodes)]
        if (length(newN) == 0)
            done <- TRUE
        else {
            oldNodes <- c(oldNodes, newN)
            numE <- length(newN)
            nedges <- AnnotationDbi::mget(newN, env = env, ifnotfound = NA)
            nedges <- nedges[!is.na(nedges)]
            oldEdges <- c(oldEdges, nedges)
            if (length(nedges) > 0)
                newN <- sort(unique(unlist(nedges)))
            else newN <- NULL
        }
    }
    rE <- vector("list", length = length(oldNodes))
    names(rE) <- oldNodes
    rE[names(oldEdges)] <- oldEdges
    rE <- lapply(rE, function(x) match(x, oldNodes))
    names(oldNodes) = oldNodes
    return(new("graphNEL", nodes = oldNodes, edgeL = lapply(rE,
        function(x) list(edges = x)), edgemode = "directed"))
}

getGOGraph<-function(term, prune=Inf){
	if(!require(graph))
		stop("Package graph is required for function getGOGraph")
	if(!exists("GOSimEnv")) initialize()	
	ontology<-get("ontology",envir=GOSimEnv)		
	if(ontology == "BP")
		G<-GOGraph(term,GOBPPARENTS)
	else if(ontology == "MF")
		G<-GOGraph(term,GOMFPARENTS)
	else if(ontology == "CC")
		G<-GOGraph(term,GOCCPARENTS)
	else
		stop(paste("ontology", ontology, "not known!"))			
	if(prune != Inf){
		dis = johnson.all.pairs.sp(G)		
		inc = unique(unlist(sapply(term, function(t) names(dis[t,])[dis[t,] < prune])))
		G = subGraph(nodes(G)[inc], G)
	}		
	G
}

calcICs<-function(DIR="."){	
	require(GO.db)
	if(!require(annotate))
		stop("Package annotate is required for function calcICs")
	if(!exists("GOSimEnv")) initialize()
	evidences<-get("evidences", envir=GOSimEnv)
	ontology<-get("ontology",envir=GOSimEnv)
	organism = get("organism", envir=GOSimEnv)
	message(paste("calculating information contents for ontology ", ontology, " using evidence codes '", paste(evidences,collapse=", "), "' (", organism,") ...",sep=""))	
	ids<- toTable(GOTERM)	 
	ids = unique(ids[ids[,"Ontology"] == ontology,"go_id"]) # these are all GO terms, which belong to the corrseponding category
	offspring <- getOffsprings()	
	gomap <- get("gomap",envir=GOSimEnv)		
	goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) # all GO terms appearing in an annotation
	goterms<-goterms[goterms %in% ids] # this is to ensure that we only get GO terms mapping to the given ontology
	tab<-table(goterms)
	na<-names(tab)			
	s<-setdiff(ids, na)  #ensure that GO terms NOT appearing in the annotation have 0 frequency
	m<-double(length(s))
	names(m)<-s
	tab<-as.vector(tab)
	names(tab)<-na
	tab<-c(tab,m)	
	ta<-sapply(ids,function(x){ t=tab[unlist(offspring[x])]; tab[x]+sum(t[!is.na(t)])})		
	names(ta)<-ids	
	IC<- -log(ta/sum(tab))	# ACHTUNG: hier muß tab und nicht ta stehen: Die Vorkommenshäufigkeit des Wurzelknotens ist die Summe der Vorkommenshäufigkeiten aller Knoten OHNE Aufsummieren der Kinder!
# # 	IC[IC == Inf] = 0 # WRONG: GO terms which are not annotated have Inf information content (NOT 0: They cannot be treated like root!!!)
	save(IC,file=file.path(DIR, paste("ICs",ontology,organism,paste(evidences,collapse="_"),".rda",sep="")))	
	message("done")			
}

getMinimumSubsumer<-function(term1, term2){	
	if(!exists("GOSimEnv")) initialize()
	ancestor<-get("ancestor",envir=GOSimEnv)
	if(term1 == term2){
		ms<-term1	
		return(ms)
	}
	an1<-unlist(ancestor[names(ancestor) == term1])
	an2<-unlist(ancestor[names(ancestor) == term2])
	case1<-which(an2 == term1)  # term1 is the ms of term2
	case2<-which(an1 == term2) # term2 is the ms of term1	  
	if(length(case1) > 0){
		ms<-term1	
	} else if(length(case2) > 0) {
		ms<-term2	
	} else {
		# find common ancestor with maximal information content
		anall<-intersect(an1, an2) 
		IC<-get("IC", envir=GOSimEnv)		
		ms<-anall[which.max(IC[anall])]
	}	
	if(is.null(ms) | length(ms) == 0)
		ms <- "NA"	
	ms
}

# get FuSSiMeg density factor
getDensityFactor<-function(term){
	if(!exists("GOSimEnv")) initialize()
	nchildren<-get("nchildren",envir=GOSimEnv)
	nparents<-get("nparents",envir=GOSimEnv)
	e<-nchildren[term] + nparents[term]	  
	betaParam<-get("betaParam",envir=GOSimEnv)
	E<-(1-betaParam)*get("Eavg",envir=GOSimEnv)/e + betaParam
	E
}

# get FuSsiMeg depth factor
getDepthFactor<-function(term,G){	
	if(!require(RBGL))
		stop("Package RBGL is required for function getDepthFactor")
	if(!exists("GOSimEnv")) initialize()	
	d<-sp.between(G,term,"all")[[1]]$length + 1  # start with depth = 1!
    D<-((d+1)/d)^get("alphaParam",envir=GOSimEnv)
    D
}

# compute FuSSiMeg enriched term similarity
getEnrichedSim<-function(term1, term2){   
	if(!require(RBGL))
		stop("Package RBGL is required for function getDepthFactor")
	if(!exists("GOSimEnv")) initialize() 
	ms<-getMinimumSubsumer(term1,term2)
	IC<-get("IC", envir=GOSimEnv)	
	if(term1 != term2){    	    	
		G<-getGOGraph(c(term1,term2))
		if(term1 != ms){
			path1=sp.between(G,term1,ms)[[1]]$path # path to ms                
			len<-length(path1)
			delta1<-sum(sapply(path1[2:len],getDepthFactor,G)*sapply(path1[2:len],getDensityFactor)*(-diff(IC[path1])))
		}
		else
			delta1<-0
		if(term2 != ms){
			path2<-sp.between(G,term2,ms)[[1]]$path # path to ms    	
			len<-length(path2)
			delta2<-sum(sapply(path2[2:len],getDepthFactor,G)*sapply(path2[2:len],getDensityFactor)*(-diff(IC[path2])))
		}
		else
			delta2<-0
		delta<-delta1 + delta2		
		sim<-1 - min(1, delta)
	}
	else
		sim<-1 
	sim<-sim * IC[term1] * IC[term2]  # correction given in equation (11) of the FuSSiMeg paper
	sim[is.na(sim)] = 0
	names(sim)<-c()   
	sim
}

# get GraSM disjunctive ancestors of a set of terms with ancestors an
getDisjAnc<-function(term, an){		
	if(!require(RBGL))
		stop("Package RBGL is required for function getDisjAnc")
	G<-getGOGraph(term)	
	disan<-matrix(0,ncol=2,nrow=0)
	for(n1 in 1:length(an)){
		if(n1 > 1){
			for(n2 in 1:(n1-1)){
				if(!separates(term, an[n1], an[n2], G) && !separates(term, an[n2], an[n1], G))
					disan<-rbind(disan,c(an[n1], an[n2]))
			}
		}
	}
	disan
}

# get GraSM common disjunctive ancestors of two terms
getDisjCommAnc<-function(term1, term2){
	if(!exists("GOSimEnv")) initialize()
	ancestor<-get("ancestor",envir=GOSimEnv)
	IC<-get("IC", envir=GOSimEnv)
	if(term1 == term2){
		return(term1)
	}
	else{
		an1<-unlist(ancestor[names(ancestor) == term1])
		an2<-unlist(ancestor[names(ancestor) == term2])
		case1<-which(an2 == term1)  # term1 is an ancestor of term2
		case2<-which(an1 == term2) # term2 is an ancestor of term1
		if(length(case1) > 0){
			ancommon<-an1
			andisj<-getDisjAnc(term1, an1)
		}
		else if(length(case2) > 0){			
			ancommon<-an2
			andisj<-getDisjAnc(term2, an2)			
		}
		else{
			ancommon<-intersect(an1,an2)
			andisj<-getDisjAnc(c(term1,term2), ancommon) # we only need to calculate the disjunctives among the common ancestors!
			#andisj = unique(rbind(getDisjAnc(term1, an1), getDisjAnc(term2, an2)))
		}		
		djca<-c()
		cond1<-sapply(ancommon, function(x) setdiff(ancommon[which(IC[ancommon] >= IC[x])],x)) # which common ancestors are more informative than a_i?
		if(length(cond1) > 0){ # look for those that are disjunctive
			names(cond1)<-ancommon					
			for(i in 1:length(cond1)){
				res<-sapply(cond1[i][[1]],function(x) any(andisj[,1] == names(cond1[i]) & andisj[,2] == x) | any(andisj[,2] == names(cond1[i]) & andisj[,1] == x))				
				if(all(res))
					djca<-c(djca, names(cond1[i]))
			}
			djca= unique(djca)			
		}	
		if(length(djca)==0)
			djca<-ancommon[which.max(IC[ancommon])]# take minimum subsumer otherwise		
		return(djca)		
	}	
}
	
# get GraSM similarity of common disjunctive ancestors of two terms
getDisjCommAncSim<-function(term1, term2, method="JiangConrath"){
	if(!exists("GOSimEnv")) initialize()
	IC<-get("IC", envir=GOSimEnv)
	djca<-getDisjCommAnc(term1, term2)	
	ICdjca<-IC[djca]	
	ICdjca<-ICdjca[!is.na(ICdjca)]							
	ICshare<-mean(ICdjca)		
	if(method == "JiangConrath")
		return(1-min(1,-2*ICshare + IC[term1] + IC[term2]))
	else if(method == "Resnik")
		return(ICshare)
	else if(method == "Lin")
		return(2*ICshare/(IC[term1]+IC[term2]))
	else
		stop(paste("getDisjCommAnc: Unknown term similarity",method))
}

# graph information content similarity related to Tanimoto-Jacard index
getGIC = function(term1, term2){
	if(!exists("GOSimEnv")) initialize()	
	if(term1 == term2){
		return(1)
	}
	IC<-get("IC", envir=GOSimEnv)
	ancestor<-get("ancestor",envir=GOSimEnv)
	an1<-unlist(ancestor[names(ancestor) == term1])
	an2<-unlist(ancestor[names(ancestor) == term2])
	ancommon = intersect(an1, an2)
	anunion = union(an1, an2)
	return(sum(IC[ancommon]) / sum(IC[anunion]))
}

# basic term similarity for index i and j in GOIDs ids
calcTermSim<-function(ids, i, j, method="JiangConrath"){	
  	calcTermSim(ids[i],ids[j], method)
}

# basic term similarity between term1 and term2
calcTermSim<-function(term1, term2, method="JiangConrath", verbose=FALSE){
	if(!exists("GOSimEnv")) initialize()
	IC<-get("IC", envir=GOSimEnv)
	if(verbose)
		message(paste("Terms:",term1,",",term2,"( method:",method,")"))	
	if(method== "Resnik")
		return(IC[getMinimumSubsumer(term1,term2)] / max(IC[IC != Inf]))   
	else if(method == "JiangConrath")
		return(1 - min(1, -2*IC[getMinimumSubsumer(term1,term2)] + IC[term1] + IC[term2]) )	
	else if(method == "Lin"){
		res = 2*IC[getMinimumSubsumer(term1,term2)]/(IC[term1]+IC[term2])
		return(ifelse(is.na(res), 0, res))
	}
	else if(method== "CoutoEnriched")
		return(getEnrichedSim(term1, term2))
	else if(method == "CoutoResnik")  
		return(getDisjCommAncSim(term1,term2, "Resnik"))
	else if(method == "CoutoJiangConrath")  
		return(getDisjCommAncSim(term1,term2, "JiangConrath"))
	else if(method == "CoutoLin"){
		res = getDisjCommAncSim(term1,term2, "Lin")
		return(ifelse(is.na(res), 0, res))
	}
	else if(method == "simIC"){ # Li et al.
		MICA = getMinimumSubsumer(term1,term2)
		res = 2*IC[MICA]/(IC[term1]+IC[term2]) * (1 - 1/(1 + IC[MICA]))
		return(ifelse(is.na(res), 0, res))
	}
	else if(method == "GIC") # graph information content
		return(getGIC(term1, term2))
	else if(method == "relevance"){ # Schlicker et al.
		MICA = getMinimumSubsumer(term1,term2)
		res = (2*IC[MICA]/(IC[term1]+IC[term2]))*(1 - exp(-IC[MICA]))
		return(ifelse(is.na(res), 0, res))
	}
	else if(method == "diffKernel"){
		K = mget("K", envir=GOSimEnv, ifnotfound=list(function(x) stop(paste("Diffusion kernel not loaded!\nPlease invoke load.diffusion.kernel().", sep = ""), call. = FALSE)))$K
		return(K[term1, term2])
	}
	else
		stop(paste("calcTermSim: Unknown term similarity",method))
}

# calculate term similarities for a list of GO terms
getTermSim<-function(termlist, method="relevance", verbose=FALSE){
	S<-matrix(0,nrow=length(termlist),ncol=length(termlist))
	colnames(S)<-termlist
	rownames(S)<-termlist
	if(method %in% c("diffKernel")){		
		K = mget("K", envir=GOSimEnv, ifnotfound=list(function(x) stop(paste("Diffusion kernel not loaded!\nPlease invoke load.diffusion.kernel().", sep = ""), call. = FALSE)))$K
		Ktmp = matrix(NA, ncol=length(termlist), nrow=length(termlist))
		diag(Ktmp) = 1
		termlist = intersect(termlist, colnames(K))
		K = K[termlist, termlist]		
		return(K)
	}
	for(i in 1:length(termlist)){
		S[i,i] <- calcTermSim(termlist[i],termlist[i], method, verbose)				
		if(i > 1){
			for(j in 1:(i-1)){				
				S[i,j]<- calcTermSim(termlist[i],termlist[j], method, verbose)				
				S[j,i]<-S[i,j]
			}
		}
	}
	S
}