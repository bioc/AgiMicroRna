\name{filterMicroRna}
\alias{filterMicroRna}
\title{ Filtering Genes }
\description{
	 Filter genes out according to their Quality Flag	
}
\usage{
filterMicroRna(ddNORM, 
	dd,
	control,
	IsGeneDetected,
	wellaboveNEG,
	limIsGeneDetected,
	limNEG,
	makePLOT=FALSE,
	targets,
	verbose=FALSE,
  writeout=FALSE)
}
\arguments{
  \item{ddNORM}{RGList with the Total Gene Signal in log2 scale to be FILTERED out according to a Quality FLAG }
  \item{dd}{RGList, containing the output from \code{readMicroRnaAFE}}
  \item{control}{logical, if \code{TRUE} it removes controls }
  \item{IsGeneDetected}{logical, if \code{TRUE} it filters genes according to gIsGeneDetected Flag. 
			Flag = 1, then gene is detected 			
 }
  \item{wellaboveNEG}{logical, if \code{TRUE} it filter genes whose expression is not above a limit
			value defined by the expression of negative controls.  
			Limit= Mean(negative) + 1.5 x sd(negative) 
 }
  \item{limIsGeneDetected}{ for a given feature xi accros samples, is the
                minimum % of features that is demanded to remain
                in at least one experimental condition with a IsGeneDetected-FLAG = 1 (Is Detected)
 }
  \item{limNEG}{for a given feature xi accros samples, is the
                minimum % of features that is demanded to remain
                in at least one experimental condition with intensity > Limit established
                for negative controls (Mean + 1.5 x SD)
 }
  \item{makePLOT}{logical, if \code{TRUE} makes QC plots with the remaining signals }
  \item{targets}{data.frame with the targets structure }
  \item{verbose}{logical, if \code{TRUE} prints out output} 
   \item{writeout}{logical, if \code{TRUE} writes out output} 
}
\details{
	Agilent Feature Extraction software provides a flag for each spot
  that identifies different quantification errors of the signal.
  Quantification flags were used to filter out signals that did not reach
  a minimum established criterion of quality. 
}
\value{
  The function returns a RGList containing the FILTERED data. In order to allow
  the tracking of those microRNAs that may have been filtered out from
  the original raw data, the following files are given:

	 NOCtrl\_exprs.txt: Log2 Normalized Total Gene Signals for the Non Control Genes 
	 NOCtrl\_FlagIsGeneDetected.txt: IsGeneDetected Flag for the Non Control Genes, 1 = detected 
	 IsNOTGeneDetected.txt: Genes that not are not detected according to IsGeneDetected Flag 

}
\references{ 
Agilent Feature Extraction Reference Guide
        \url{http://www.Agilent.com} 
}
\author{ Pedro Lopez-Romero }
\examples{
data(dd.micro,verbose=FALSE)
data(targets.micro,verbose=FALSE)
ddTGS=tgsMicroRna(dd.micro,half=TRUE,makePLOT=FALSE,verbose=FALSE)
ddNORM=tgsNormalization(ddTGS,'quantile',
                      makePLOTpre=FALSE,makePLOTpost=TRUE,targets.micro,verbose=FALSE)
ddPROC=filterMicroRna(ddNORM,
		      dd.micro,
                      control=TRUE,
                      IsGeneDetected=TRUE,
                      wellaboveNEG=FALSE,
                      limIsGeneDetected=50,
                      limNEG=25,
                      makePLOT=FALSE,
                      targets.micro,
		      verbose=TRUE)
}
\keyword{documentation}
\keyword{utilities}

