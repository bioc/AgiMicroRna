\name{tgsMicroRna}
\alias{tgsMicroRna}
\title{ Getting the Total Gene Signal }
\description{
	The function creates an uRNAList containing the TotalGeneSignal computed
	by the Agilent Feature Extraction software. This signal can be 
	used for the statistical analysis after a possible normalization step.  
}
\usage{
tgsMicroRna(dd, offset, half, makePLOT=FALSE, verbose=FALSE)
}
\arguments{
  \item{dd}{uRNAList, containing the output from \code{readMicroRnaAFE} }
  \item{offset}{integer. To use this option set half = FALSE }
  \item{half}{logical, if \code{TRUE} half option is used }
  \item{makePLOT}{logical, if \code{TRUE} QC plots with the Total Gene Signal are displayed  }
  \item{verbose}{logical, if \code{TRUE} prints out some summary results }
}
\details{
	The function creates a uRNAList oobject that contains in the uRNAList$TGS, uRNAList$TPS, 
	uRNAList$meanS & uRNAList$procS the Total Gene Signal (TGS) 
	as computed by the Agilent Feature Extraction algorithms. This TGS is not in
	log2 scale. All the replicated genes have
	the same estimated TGS, and the function simply picks one gene from each set 
	of replicated genes. To mantain the format of the uRNAList, every selected gene 
	retains a probe name attach to them. This probe name is not meaningful any more,  
	since the signal corresponds to the total gene signal and not to the probe signal. 
	The TGS processed by AFE contains some negative values. To get signals with positive values 
	we can either add a positive small constant to all the signals (offset) or	
	we can select the 'half' option, which set to 0.5 all the values that are smaller than
	0.5. To use the offset option we have to set half=FALSE, otherwise the half method
	is used by default. The offset option, adds to each signal the quantity (abs( min(ddTGS$TGS)) + offset),
	where ddTGS$TGS is the matrix that contains the TotalGeneSignal. 
}
\value{
	uRNAList containing the TotalGeneSignal computed
	by the Agilent Feature Extraction software 
	Optionally, it can generate a boxplot, a density plot and a MA plot
	with the Total Gene Signal. 
}
\references{ 
	 Agilent Feature Extraction Reference Guide
         \url{http://www.Agilent.com}
 }
\author{ Pedro Lopez-Romero }
\examples{
	data(dd.micro)
	data(targets.micro)
	ddTGS=tgsMicroRna(dd.micro,half=TRUE,makePLOT=FALSE,verbose=FALSE)
}
\keyword{documentation}
\keyword{utilities}

