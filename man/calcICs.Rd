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
   new version of the "GO" library has been installed.
   By default the information contents are precomputed using all
   evidence codes and evidence codes "IMP, IGI, IDA, IEP, IPI" together.   
}
\value{
  Puts a file named "ICs<ontology><evidence levels>.rda" in the current
  directoy. This file has then to be moved into the data directory of
  gosim. 
}
\seealso{\code{\link{setEvidenceLevel}}}
\examples{		
\dontrun{
  setOntology("CC")
  calcICs()
} # --> this may take some time ...
}
\keyword{ file }% at least one, from doc/KEYWORDS
