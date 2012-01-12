\name{calcICs}
\alias{calcICs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate information contents of GO terms.}
\description{
   Recalculates the information content of all GO terms.
}
\usage{
calcICs(DIR=".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{DIR}{ directory where to put the resulting files}
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
  Puts a file named "ICs<ontology><organism><evidence levels>.rda" in 
  directoy DIR. It can be used afterwards by calling \code{\link{setOntology}}.   
}
\seealso{\code{\link{setEvidenceLevel}}}
\examples{		
\dontrun{  	
	setEvidenceLevel("IMP")
	setOntology("CC", loadIC=FALSE) # important: setOntology assumes that the IC file already exists. To prevent an error message we need the second argument
	calcICs()  
} # --> this may take some time ...
}
\keyword{ file }% at least one, from doc/KEYWORDS
