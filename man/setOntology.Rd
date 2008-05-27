\name{setOntology}
\alias{setOntology}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set an ontology as base for subsequent computations.}
\description{
   Sets the ontology that all subsequent computations are based on and
   loads the information content of all GO terms within this
   ontology. At load time of the library the default ontology is
   "BP". Furtherm, on running this function the environment GOSimEnv is
   reinitialized, i.e. all global settings or parameters used in the
   library are reset to their default values. 
}
\usage{
setOntology(ont = "BP")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ont}{ the ontology to use ("BP","MF","CC") }
}
\details{
  The following ontologies can be used:
\item{"BP"}{ biological process}
\item{"MF"}{ molecular function}
\item{"CC"}{ cellular component}
}
\value{
 none.
}
\author{ Holger Froehlich }
\examples{
 # set ontology to "molecular function"
\dontrun{
 setOntology("MF")
 # calculate Resnik similarity of two GO terms within this ontology
 getTermSim(c("GO:0004060","GO:0003867"),method="Resnik")
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
