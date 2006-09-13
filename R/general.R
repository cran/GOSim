initialize<-function(){
	print("initializing GOSim package ...")
# 	require(GOstats)  
#   	require(mclust) 
#   	require(cluster) 
	assign("GOSimEnv",new.env(parent=globalenv()),envir=.GlobalEnv)  
  	setEvidenceLevel("all")
  	setOntology("BP")
  	print("finished.")
}

getOffsprings<-function(){
  if(!exists("GOSimEnv")) initialize()
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-as.list(GOBPOFFSPRING)
  else if(ontology == "MF")
    res<-as.list(GOMFOFFSPRING)
  else if(ontology == "CC")
    res<-as.list(GOCCOFFSPRING)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

getAncestors<-function(){
  if(!exists("GOSimEnv")) initialize()
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-as.list(GOBPANCESTOR)
  else if(ontology == "MF")
    res<-as.list(GOMFANCESTOR)
  else if(ontology == "CC")
    res<-as.list(GOCCANCESTOR)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

getParents<-function(){
  if(!exists("GOSimEnv")) initialize()
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-as.list(GOBPPARENTS)
  else if(ontology == "MF")
    res<-as.list(GOMFPARENTS)
  else if(ontology == "CC")
    res<-as.list(GOCCPARENTS)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

getChildren<-function(){
  if(!exists("GOSimEnv")) initialize()
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-as.list(GOBPCHILDREN)
  else if(ontology == "MF")
    res<-as.list(GOMFCHILDREN)
  else if(ontology == "CC")
    res<-as.list(GOCCCHILDREN)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

# filter GO mapping for given evidence levels
setEvidenceLevel<-function(evidences="all"){	
        if(!exists("GOSimEnv")) initialize()	
	assign("evidences", evidences, envir=GOSimEnv)	
	gomap<-as.list(GOENTREZID2GO)	
	if((length(evidences) > 1) || (evidences!="all")){
		gomap<-sapply(gomap,function(x) sapply(x,function(y) y$Evidence %in% evidences))
		gomap<-sapply(gomap, function(x) x[which(x)])
		gomap<-gomap[sapply(gomap,function(x) length(x) > 0)]
	}
	assign("gomap", gomap, envir=GOSimEnv)
}

setOntology<-function(ont="BP"){
	if(!exists("GOSimEnv")) initialize()
	assign("ontology", ont, envir=GOSimEnv)		
	ontology<-get("ontology",envir=GOSimEnv)
	evidences<-get("evidences",envir=GOSimEnv)	
	data(list=paste("ICs",ontology,evidences,sep=""),envir=GOSimEnv)	
	IC<-get("IC",envir=GOSimEnv)
	IC<-IC/max(IC[IC!=Inf])
	IC["all"]=0
	assign("IC", IC, envir=GOSimEnv)
 	assign("ancestor", getAncestors(), envir=GOSimEnv) 	
 	assign("children", getChildren(), envir=GOSimEnv)
 	assign("parents", getParents(), envir=GOSimEnv) 
 	children<-get("children",envir=GOSimEnv) 	
 	parents<-get("parents",envir=GOSimEnv) 	
 	assign("nchildren", sapply(children,length) , envir=GOSimEnv)
 	assign("nparents", sapply(parents,length), envir=GOSimEnv)
 	nchildren<-get("nchildren",envir=GOSimEnv) 	
 	nparents<-get("nparents",envir=GOSimEnv) 	
 	assign("Eavg", (sum(nchildren) + sum(nparents))/ length(union(names(nchildren), names(nparents))), envir<-GOSimEnv)
  	assign("alphaParam", 0.5, envir=GOSimEnv)
  	assign("betaParam", 0.5, envir=GOSimEnv)  	    	
}