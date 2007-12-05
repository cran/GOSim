# perform a GO gene set enrichment analysis for a specific cluster using topGO
GOenrichment = function(allpvalues, fdr=0.05, cutoff=0.01){	
	require(topGO)
	require(annotate)
	ontology = get("ontology", envir=GOSimEnv)
	gomap <- get("gomap",env=GOSimEnv)	
	anno <- gomap[names(allpvalues)]
	goterms <- sapply(anno,function(x) names(x))
	topdiffgenes <- function(allScore) {
       		return(allScore < fdr)
     	}
	GOdata = new("topGOdata", ontology = ontology, allGenes = allpvalues, geneSel = topdiffgenes, annot = annFUN.gene2GO, gene2GO = goterms) 
	test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "KS test", cutOff = 0.01) 
	res = getSigGroups(GOdata, test.stat) 
	sigterms = score(res)[score(res) < cutoff]
	goids<-as.list(GOTERM)	
	goids<-goids[sapply(goids, function(x) Ontology(x) == ontology)]
	list(GOTerms=goids[names(sigterms)], p.values = sigterms)
}