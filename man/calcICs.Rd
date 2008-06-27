\name{calcICs}
\alias{calcICs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate information contents of GO terms.}
\description{
   Recalculates the information content of all GO terms.
}
\usage{
calcICs()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  none
}
\details{
   This functions should only be invoked, if one wants to calculate the
   information content for GO terms with respect to combinations of
   evidence codes other than the precomputed ones or, if a 
   new version of the organism annotation packages has been installed.
   By default the information contents are precomputed using all
   evidence codes and evidence codes "IMP, IGI, IDA, IEP, IPI" together.   
}
\value{
  Puts a file named "ICs<ontology><evidence levels>.rda" in the current
  directoy. This file has then to be moved into the data directory of
  GOSim. It can be used afterwards by calling \code{\link{setOntology}}. 
  Sometimes it is necessary to switch to the GOSim directory for this purpose.
}
\seealso{\code{\link{setEvidenceLevel}}}
\examples{		
\dontrun{
  if(require("annotate")){	
	setEvidenceLevel("all")
	setOntology("CC") # important: setOntology assumes that the IC file already exists. To prevent an error message we need this workaround
	setEvidenceLevel("IMP")
	calcICs()
  }
} # --> this may take some time ...
}
\keyword{ file }% at least one, from doc/KEYWORDS
