\name{cvArray}
\alias{cvArray}
\title{ Coefficient of variation of replicated probes within array}
\description{
	Identifies replicated features at probe and at gene level and 
	computes the coefficient of variation of the array
}
\usage{
cvArray(ddDUP, foreground = c("MeanSignal", "ProcessedSignal"),targets,verbose=FALSE)
}
\arguments{
  \item{ddDUP}{RGList, containing the output from \code{readMicroRnaAFE}  }
  \item{foreground}{Specifies the signal used, only "MeanSignal" or "ProcessedSignal" can be used}
  \item{targets}{data.frame with the target structure  }
  \item{verbose}{logical, if \code{TRUE} prints out output}	
}
\details{
  In the Agilent microRNA platforms the features are replicated 
  at a probe level and normally, a single microRNA is interrogated
  by either two or four sets of replicated probes. The replication of the probes allows 
  computing the coefficient of variation (CV) for each array as a measure of 
  the reproducibility of the array. The CV is computed for every set of 
  replicated probes and the CV median is reported as the array CV. A lower array CV 
  indicates a better array reproducibility.
}
\value{
	It prints out the results of the replication for the NON CONTROL FEATURES 
	at a probe and gene level. 
}
\author{ Pedro Lopez-Romero }
\examples{
\dontrun{
data(dd.micro)
data(targets.micro)

cvArray(dd.micro,"MeanSignal",targets.micro,verbose=TRUE)

graphics.off()
}
}
\keyword{documentation}
\keyword{utilities}
