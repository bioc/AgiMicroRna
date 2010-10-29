\name{mvaMicroRna}
\alias{mvaMicroRna}
\title{MA plot}
\description{
	For each array, the M value is computed for every spot as the difference
        between the spot intensity in the array and the median intensity for that feature over the whole
        set of arrays.
        Every kind of feature is identified with different color (microRNA genes, positive controls, etc ...)
	The input must be an uRNAList object created by the user, in such a way that the uRNAList$meanS field 
	contains the expression matrix that we want to use in log2 scale (see example below)  
	The gProcessedSignal computed by the Agilent Feature Extaction software normally
	contains negative values, so a small constant has to be added to the signals before log 
	tranformation. 
}
\usage{
mvaMicroRna(uRNAList, maintitle, verbose=FALSE)
}
\arguments{
  \item{uRNAList}{A uRNAList object. It uses the expression matrix stored in the uRNAList$meanS slot.
		Input expression matrix should be in log2 scale}
  \item{maintitle}{ character to indicate the title of the graph }
  \item{verbose}{ logical, if \code{TRUE}it prints details }
}
\author{ Pedro Lopez-Romero }
\examples{
data(dd.micro)
op=par(mfrow=c(1,1),ask=TRUE)

MMM=dd.micro$procS  		  ## gProcessedSignal 

min=min(MMM)        		  ## transforming gProcessedSignal to positive values 
for(i in 1:dim(MMM)[2]){    	  ## before log2 transformation 		
MMM[,i]=MMM[,i]+(abs(min)+ 5)	
}
        ddaux=dd.micro
        ddaux$meanS=log2(MMM)        
	mvaMicroRna(ddaux,maintitle="ProcessedSignal",verbose=FALSE)
	rm(ddaux)
par(op)
}
\keyword{documentation}
\keyword{utilities}

