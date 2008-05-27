\name{setEvidenceLevel}
\alias{setEvidenceLevel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specifies to use only GO terms with given evidence codes.}
\description{
   Specifies to use only GO terms with given evidence codes. 
}
\usage{
setEvidenceLevel(evidences = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evidences}{ character vector of evidence codes }
}
\details{
   Each evidence code can be one of:
\item{"IMP"}{inferred from mutant phenotype}
\item{"IGI"}{inferred from genetic interaction}
\item{"IPI"}{inferred from physical interaction}
\item{"ISS"}{inferred from sequence similarity}
\item{"IDA"}{inferred from direct assay}
\item{"IEP"}{inferred from expression pattern}
\item{"IEA"}{inferred from electronic annotation}
\item{"TAS"}{traceable author statement}
\item{"NAS"}{non-traceable author statement}
\item{"ND"}{no biological data available}
\item{"IC"}{inferred by curator}
Entrez Gene ids for which no GO associations exist are left out of the environment.

Mappings were based on data provided by:

Entrez Gene:
<URL: http://gopher5/compbio/annotationSourceData/ftp.ncbi.nlm.nih.gov/gene/DATA/>.
Built: Source data downloaded from Entrez Gene on Fri Sep 30 02:51:32 2005

}
\value{
  none.
}
\references{<URL: http://www.ncbi.nih.gov/entrez/query.fcgi?db=gene>}
\author{ Holger Froehlich }
\note{ The method directly uses the "GO" library to obtain the mappings
  of the Entrez Gene IDs to GO terms. Internally this mapping is
  especially used for calculating the information content of the GO
  terms. By default all evidence codes are used for this. If another
  behavior is wanted, one has to recalculate the information content of
  all GO terms via \code{\link{calcICs}}.
  The evidence level also influences the behavior of
  \code{\link{filterGO}} and \code{\link{getGOInfo}}. 
}
\seealso{ \code{\link{calcICs}}, \code{\link{filterGO}}, \code{\link{getGOInfo}} }
\examples{
	\dontrun{setEvidenceLevel("all")} # the default behavior	
}
\keyword{ file }% at least one, from doc/KEYWORDS

