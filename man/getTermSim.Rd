\name{getTermSim}
\alias{getTermSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get pairwise GO term similarities.  }
\description{
   Returns the pairwise similarities between GO terms. Different
   calculation method are implemented.
}
\usage{
getTermSim(termlist, method = "relevance", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{termlist}{ character vector of GO terms }  
  \item{method}{ one of the supported methods for GO term similarity (see below) }
  \item{verbose}{ print out various information or not }
}
\details{
  Currently the following methods for computing GO term similarities are implemented: 
  \describe{
\item{"Resnik"}{information content of minimum subsumer (ICms) [1]}
\item{"JiangConrath"}{\eqn{1 - \min(1, IC(term1) - 2ICms + IC(term2))} [2]}
\item{"Lin"}{\eqn{\frac{2ICms}{(IC(term1) + IC(term2))}} [3]}
\item{"CoutoEnriched"}{FuSSiMeg enriched term similarity by Couto et al. [4]. Requires enrichement factors to be set by \code{\link{setEnrichmentFactors}}}.
\item{"CoutoResnik"}{average information content of common disjunctive ancestors of \code{term1} and \code{term2} (ICshare) [5]}
\item{"CoutoJiangConrath"}{\eqn{1 - \min(1, IC(term1) -2ICshare + IC(term2))} [5]}
\item{"CoutoLin"}{\eqn{\frac{2ICshare}{(IC(term1) + IC(term2))}} [5]}
\item{"diffKernel"}{diffusion kernel similarity from a pre-loaded kernel matrix (see \code{\link{load.diffusion.kernel}}). The diffusion kernel is calculated using one of the methods described in [6].}
\item{"relevance"}{sim_Lin * (1 - exp(-ICms)) [7]}
\item{"GIC"}{summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 [8]}
%%\item{simIC}{sim_Lin * (1 - 1/(1 + ICms)) [9]}
}
}
\value{
n x n matrix (n = number of GO terms) with similarities between GO terms scaled to [0,1]. If a GO term does not exist for the currently set ontology, the similarity is set to "NA".
}
\references{
[1] P. Resnik, Using Information Content to evaluate semantic similarity in a taxonomy, Proc. 14th Int. Conf. Artificial Intel., 1995 

[2] J. Jiang, D. Conrath, Semantic Similarity based on Corpus Statistics and Lexical Taxonomy, Proc. Int. Conf. Research in Comp. Ling., 1998 

[3] D. Lin, An Information-Theoretic Definition of Similarity, Proc. 15th Int. Conf. Machine Learning, 1998 

[4] F. Couto, M. Silva, P. Coutinho, Implementation of a Functional Semantic Similarity Measure between Gene-Products, DI/FCUL TR 03-29, Department of Informatics, University of Lisbon, 2003 

[5] Couto, F.; Silva, M. & Coutinho, P., Semantic Similarity over the Gene Ontology: Family Correlation and Selecting Disjunctive Ancestors, Conference in Information and Knowledge Management, 2005  

[6] Lerman G. & Shaknovich B., Defining Functional Distance using Manifold Embeddings of Gene Ontology Annotations, PNAS, 104(27): 11334 - 11339, 2007

[7] A. Schlicker, F. Domingues, J. Rahnenfuehrer, T. Lengauer, A new measure for functional similarity of gene products based on Gene Ontology, BMC Bioinformatics, 7, 302, 2006.

[8] C. Pesquita, D. Faria, H. Bastos, A. Falcao, F. Couto, Evaluating GO-based Semantic Similarity Measures, In: Proc. 10th Annual Bio-Ontologies Meeting 2007, 37 - 40, 2007

%%[9] B. Li, J. Wang, A. Feltus, J. Zhou, F. Luo, Effectively Integrating Information Content and Structural Relationship to Improve the GO-based Similarity Measure Between Proteins, BMC Bioinformatics, 2009, in press. 
}
\author{ Holger Froehlich }
\note{ All calculations use normalized information contents for each GO
  term. Normalization is achieved by dividing each information content
  by the maximum information content within the currently set ontology
  ("BP","MF","CC") 
}
\seealso{ \code{\link{getMinimumSubsumer}},
  \code{\link{getDisjCommAnc}}, \code{\link{setEnrichmentFactors}},
  \code{\link{setOntology}},  \code{\link{load.diffusion.kernel}} } 
\examples{
 # setOntology("BP")

 # Lin's method
 getTermSim(c("GO:0006955","GO:0007584"),method="Lin")
 if(require(graph)){
 	# Couto's method combined with Jiang-Conrath distance
 	getTermSim(c("GO:0006955","GO:0007584"),method="CoutoJiangConrath") 

 	# set enrichment factors
 	setEnrichmentFactors(alpha=0.1,beta=0.5) 
 	getTermSim(c("GO:0006955","GO:0007584"),method="CoutoEnriched")
 }
}
\keyword{ file }% at least one, from doc/KEYWORDS

