\name{calc.diffusion.kernel}
\alias{calc.diffusion.kernel}
\alias{load.diffusion.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculation and loading of diffusion kernel matrices}
\description{
  Manifold embeddings of gene ontology terms via diffusion kernel techniques. Diffusion kernels are positive semidefinite similarity measures calculated from the graph Laplacian. They are interpreted as the result of a local heat diffusion process along the graph structure. 
}
\usage{
calc.diffusion.kernel(method="diffKernelLapl", m=7)

load.diffusion.kernel(method="diffKernelLapl")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{one of "diffKernelLapl", "diffKernelpower", "diffKernelLLE", "diffKernelexpm"}
  \item{m}{(1) Half the power of the transition probability matrix (an integer > 0). (2) an arbitrary positive time constant for the exponential diffusion kernel}
}
\details{
  \describe{
\item{"diffKernelLapl"}{pseudo inverse of the (unnormalized) graph Laplacian: Takes into account all powers of diffusion and incorporates all paths from one node to another one.}
\item{"diffKernelpower"}{even power of the transition probability matrix: Takes into account local transitions of path length m}
\item{"diffKernelLLE"}{local linear embedding into an Euclidean space: The focus is to preserve local distances to nearest neighbors. The LLE kernel emphasizes short-range interactions between GO terms.}
\item{diffKernelexpm"}{expm(-mL), where t is a positive constant, L is the (unnormalized) graph Laplacian and expm denotes the matrix exponential. This kernel takes into account all positive integer powers of diffusion, but with an exponential decay of the influence of long-range interactions.}
}
}
\value{
  \code{calc.diffusion.kernel} puts a kernel matrix / similarity matrix named "<method><ontology><organism><evidence levels>.rda" in the current directoy. This file has then to be moved into the data directory of
  GOSim. It can be used afterwards by calling \code{\link{load.diffusion.kernel}}. 
  Sometimes it is necessary to switch to the GOSim directory for this purpose.
}
\seealso{\code{\link{load.diffusion.kernel}}}
\references{
Lerman G. & Shaknovich B., Defining Functional Distance using Manifold Embeddings of Gene Ontology Annotations, PNAS, 104(27): 11334 - 11339, 2007
}
\examples{		
\dontrun{
	calc.diffusion.kernel(method="diffKernelpower", m=5)
	load.diffusion.kernel("diffKernelpower")
	getTermSim(c("GO:0007166","GO:0007267","GO:0007584","GO:0007165","GO:0007186"),method="diffKernel",verbose=FALSE)
} # --> this may take some time ...
}
\keyword{ file }% at least one, from doc/KEYWORDS
