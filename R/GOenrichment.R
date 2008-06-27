# perform a GO gene set enrichment analysis using topGO
GOenrichment = function(allpvalues, fdr=0.05, cutoff=0.01){	
	if(!require(topGO) | !require(annotate))
		stop("Packages topGO and annotate required for function GOenrichment")
	ontology = get("ontology", envir=GOSimEnv)
	gomap <- get("gomap",env=GOSimEnv)	
	anno <- gomap[names(allpvalues)]
	goterms = sapply(anno, function(x) sapply(x, function(y) y$Ontology == ontology))
	goterms <-sapply(goterms, function(x) names(x[which(x)]))
	goterms <- goterms[sapply(goterms,length) > 0]
	topdiffgenes <- function(allScore) {
       		return(allScore < fdr)
     	}
	GOdata = new("topGOdata", ontology = ontology, allGenes = allpvalues, geneSel = topdiffgenes, annot = annFUN.gene2GO, gene2GO = goterms) 
	test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "KS test", cutOff = cutoff) 
	res = getSigGroups(GOdata, test.stat) 
	sigterms = score(res)[score(res) < cutoff]
	sigGOs = names(sigterms)
	genes = sapply(1:length(sigGOs), function(s) genesInTerm(GOdata, whichGO=sigGOs[s]))
	goids<- toTable(GOTERM)	
	anno = unique(goids[goids[,"go_id"] %in% names(sigterms), c("go_id","Term","Ontology","Definition")])
	list(GOTerms=anno, p.values = sigterms, genes = genes)
}
