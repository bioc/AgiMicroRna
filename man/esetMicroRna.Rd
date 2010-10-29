\name{esetMicroRna}
\alias{esetMicroRna}
\title{ ExpressionSet object from a uRNAList }
\description{
	It creates an 'ExpressionSet' object from a 'uRNAList' with unique probe names  
	Tipically, the 'uRNAList object' contains the Total Gene Processed data
}
\usage{
esetMicroRna(uRNAList, targets, makePLOT=FALSE,verbose=FALSE)
}
\arguments{
  \item{uRNAList}{An uRNAList containing normally the processed data }
  \item{targets}{data.frame with the targets structure }
  \item{makePLOT}{logical, if \code{TRUE} it makes a 'heatmap' with the 100 greater variance
                genes, a 'hierarchical cluster' with all the genes and a pca plot }
  \item{verbose}{logical, if \code{TRUE} prints out output}
}
\details{
	It creates an ExpressionSet object from a uRNAList.
        Usually this function is applied to an uRNAList object containing
        the Total Gene Processed data.
}
\value{
	 An ExpressionSet object
}
\author{ Pedro Lopez-Romero }
\seealso{
	An 'uRNAList' example containing proccesed data is in 
        \code{ddPROC} and an overview of how the processed data is produced
        is given in \code{filterMicroRna}
}
\keyword{documentation}
\keyword{utilities}
