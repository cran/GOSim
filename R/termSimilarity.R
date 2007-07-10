setEnrichmentFactors<-function(alpha=0.5, beta=0.5){
	if(!exists("GOSimEnv")) initialize()
	assign("alphaParam",alpha, envir=GOSimEnv)
	assign("betaParam",beta, envir=GOSimEnv)
}

getGOGraph<-function(term){
	if(!exists("GOSimEnv")) initialize()
	ontology<-get("ontology",env=GOSimEnv)
	require(GOstats)
	if("package:GO.db"%in%search())
		detach(package:GO.db)
	if(ontology == "BP")
		G<-GOGraph(term,GOBPPARENTS)
	else if(ontology == "MF")
		G<-GOGraph(term,GOMFPARENTS)
	else if(ontology == "CC")
		G<-GOGraph(term,GOCCPARENTS)
	else
		stop(paste("ontology", ontology, "not known!"))		
	G
}

calcICs<-function(){	
	require(GO)
	if(!exists("GOSimEnv")) initialize()
	evidences<-get("evidences", envir=GOSimEnv)
	ontology<-get("ontology",envir=GOSimEnv)
	print(paste("calculating information contents for ontology", ontology, "using evidence codes", paste(evidences,collapse=", "), "..."))	
	ids<-as.list(GOTERM)
	require(annotate)
	ids<-names(ids[sapply(ids, function(x) Ontology(x) == ontology)])	
	offspring<-getOffsprings()	
	gomap<-get("gomap",env=GOSimEnv)	
	goterms<-unlist(sapply(gomap, function(x) names(x))) # all GO terms appearing in an annotation
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
	IC<- -log(ta/sum(ta))	
	save(IC,file=paste("ICs",ontology,paste(evidences,collapse="_"),".rda",sep=""))
	print("done")			
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
	if(is.null(ms))
		ms <- "NA"
	ms
}

# get FuSSiMeg density factor
getDensityFactor<-function(term){
	if(!exists("GOSimEnv")) initialize()
	nchildren<-get("nchildren",env=GOSimEnv)
	nparents<-get("nparents",env=GOSimEnv)
	e<-nchildren[term] + nparents[term]	  
	betaParam<-get("betaParam",envir=GOSimEnv)
	E<-(1-betaParam)*get("Eavg",env=GOSimEnv)/e + betaParam
	E
}

# get FuSsiMeg depth factor
getDepthFactor<-function(term,G){	
	if(!exists("GOSimEnv")) initialize()	
	d<-sp.between(G,term,"all")[[1]]$length  	
    	D<-((d+1)/d)^get("alphaParam",envir=GOSimEnv)
    	D
}

# compute FuSSiMeg enriched term similarity
getEnrichedSim<-function(term1, term2){   
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
    names(sim)<-c()   
    sim
}

# get GraSM disjunctive ancestors of a set of terms with ancestors an
getDisjAnc<-function(term, an){		
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
# 			andisj<-unique(rbind(getDisjAnc(term1,an1), getDisjAnc(term2,an2)))			
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

# basic term similarity for index i and j in GOIDs ids
calcTermSim<-function(ids, i, j, method="JiangConrath"){	
  	calcTermSim(ids[i],ids[j], method)
}

# basic term similarity between term1 and term2
calcTermSim<-function(term1, term2, method="JiangConrath", verbose=TRUE){
	if(!exists("GOSimEnv")) initialize()
	IC<-get("IC", envir=GOSimEnv)
	if(verbose)
		print(paste("Terms:",term1,",",term2,"( method:",method,")"))	
	if(method== "Resnik")
		return(IC[getMinimumSubsumer(term1,term2)])   
	else if(method == "JiangConrath")
		return(1 - min(1, -2*IC[getMinimumSubsumer(term1,term2)] + IC[term1] + IC[term2]) )	
	else if(method == "Lin")
		return(2*IC[getMinimumSubsumer(term1,term2)]/(IC[term1]+IC[term2]))	
	else if(method== "CoutoEnriched")
		return(getEnrichedSim(term1, term2))
	else if(method == "CoutoResnik")  
		return(getDisjCommAncSim(term1,term2, "Resnik"))
	else if(method == "CoutoJiangConrath")  
		return(getDisjCommAncSim(term1,term2, "JiangConrath"))
	else if(method == "CoutoLin")  
		return(getDisjCommAncSim(term1,term2, "Lin"))
	else
		stop(paste("calcTermSim: Unknown term similarity",method))
}

# calculate term similarities for a list of GO terms
getTermSim<-function(termlist, method="JiangConrath", verbose=TRUE){
	S<-matrix(0,nrow=length(termlist),ncol=length(termlist))
	colnames(S)<-termlist
	rownames(S)<-termlist
	for(i in 1:length(termlist)){
		S[i,i]<-1
		if(i > 1){
			for(j in 1:(i-1)){				
				S[i,j]<- calcTermSim(termlist[i],termlist[j], method, verbose)				
				S[j,i]<-S[i,j]
			}
		}
	}
	S
}