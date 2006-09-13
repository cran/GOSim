\name{getGOInfo}
\alias{getGOInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Obtain GO terms and their description for a list of genes. }
\description{
  Oobtain the GO terms and their description for a list of genes. 
}
\usage{
getGOInfo(geneIDs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneIDs}{ character vector of Entrez gene IDs }
}
\value{
 List with entry names equal to the gene IDs. Each list contains a
 sublist with entry names equal to the GO terms associated to the
 corresponding gene ID. Each entry also contains a description of the GO
 term, its definition and the ontology ("BP","CC","MF") it belongs to. 
}
\author{ Holger Froehlich }
\note{ The corresponding information is directly extracted from the "GO"
  library. The result depends on the currently set ontology
  ("BP","MF","CC"), i.e. only GO terms within the actual ontology are
  considered. 
}
\seealso{ \code{\link{setOntology}} }
\examples{
 setOntology("BP")
 getGOInfo(c("207","7494"))
}
\keyword{ file }% at least one, from doc/KEYWORDS

