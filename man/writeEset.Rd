\name{writeEset}
\alias{writeEset}
\title{ Writes the expression data matrix of an ExpressionSet object in a txt file }
\description{
        Writes the expression data matrix of an ExpressionSet object in a file.
}
\usage{
writeEset(eset, ddPROC, targets,verbose=FALSE)
}
\arguments{
  \item{eset}{ An Expression object, normally containing the procesed data}
  \item{ddPROC}{An RGList object, normally containing the processed data  }
  \item{targets}{data.frame with the targets structure }
  \item{verbose}{logical, if \code{TRUE} prints out output}
}
\details{
        Writes the expression data matrix of an ExpressionSet object in a file.
}
\author{ Pedro Lopez-Romero }
\seealso{
        An 'RGList' example containing proccesed data is in
        \code{ddPROC} and an overview of how the processed data is produced
        is given in \code{filterMicroRna}. The ExpressionSet object can be
	generated using \code{esetMicroRna}
}

\examples{
\dontrun{
data(ddPROC)
data(targets.micro)
esetPROC=esetMicroRna(ddPROC,targets.micro,makePLOT=TRUE,verbose=FALSE)
writeEset(esetPROC,ddPROC,targets.micro,verbose=TRUE)
}
}
\keyword{documentation}
\keyword{utilities}

