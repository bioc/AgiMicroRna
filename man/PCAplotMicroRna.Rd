\name{PCAplotMicroRna}
\alias{PCAplotMicroRna}
\title{ PCA plot }
\description{
	It is a wrapper for the 'plotPCA' of the 'affycoretools' package
}
\usage{
PCAplotMicroRna(eset, targets)
}
\arguments{
  \item{eset}{ An Expression Set object }
  \item{targets}{data.frame with the target structure  }
}
\author{ Pedro Lopez-Romero }
\examples{
	data(targets.micro)
	data(ddPROC)
	esetPROC=esetMicroRna(ddPROC,targets.micro,makePLOT=FALSE,verbose=FALSE)

	PCAplotMicroRna(esetPROC,targets.micro)
}
\keyword{documentation}
\keyword{utilities}

