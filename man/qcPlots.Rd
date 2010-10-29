\name{qcPlots}
\alias{qcPlots}
\title{ Plots for Quality Assessment}
\description{
  It creates BoxPlots, Density Plots, MA plots, RLE plots and hierachical 
  clustering plots with the sample data set. 
}

\usage{
qcPlots(dd, 
	offset, 
	MeanSignal=TRUE, 
	ProcessedSignal=FALSE, 
	TotalProbeSignal=FALSE, 
	TotalGeneSignal=FALSE, 
	BGMedianSignal=FALSE, 
	BGUsed=FALSE, 
	targets)
}
\arguments{
  \item{dd}{A uRNAList object containing the ouput from  \code{readMicroRnaAFE}} 
  \item{offset}{numeric value to add to the intensities before log transforming}
  \item{MeanSignal}{logical, if \code{TRUE} "gMeanSignal" is used}
  \item{ProcessedSignal}{logical, if \code{TRUE} "gProcessedSignal"  is used }
  \item{TotalProbeSignal}{logical, if \code{TRUE} "gTotalProbeSignal" is used}
  \item{TotalGeneSignal}{logical, if \code{TRUE} "gTotalGeneSignal" is used}
  \item{BGMedianSignal}{logical, if \code{TRUE} "gBGMedianSignal" is used  }
  \item{BGUsed}{logical, if \code{TRUE}  "gBGUsed" is used }
  \item{targets}{data.frame with the target structure  }
}
\details{
  The signals loaded from the AFE data files can be used for the quality 
  assesment using the graphical utilities included in the \code{qcPlots} function. 
  For the gMeanSignal, the  BoxPlots, Density Plots, MA plots, RLE plots
  and hierachical clustering plots are done. For the gProcessedSignal the same 
  plots are done, except the hierarchical clustering. For the gTotalProbeSignal 
  and the gTotalGeneSignal only the  BoxPlots and Density Plots are done, and 
  finally, for the Background signals only the Boxplots are done.   
	
}
\references{ 
	Boldstad B.M., Collin F., Brettschneider J., Simpson, K.,
        Cope L., Irizarry R. A., Speed T. P. Quality Assesement of
        Affymetrix GeneChip Data. In Bioinformatics and Computational
        Biology Solutions Using R and Bioconductor. (eds.)  Gentleman R.,
        Carey V. J., Huber W., Irizarry R. A., Dudoit S. (2005). Springer.
 }
\author{ Pedro Lopez-Romero }
\seealso{ \code{boxplotMicroRna},\code{plotDensityMicroRna},\code{RleMicroRna}, \code{mvaMicroRna} 
	and \code{hierclusMicroRna} }
\examples{
\dontrun{
data(dd.micro)
qcPlots(dd.micro,offset=5,
       		MeanSignal=TRUE,
        	ProcessedSignal=TRUE,
        	TotalProbeSignal=TRUE,
        	TotalGeneSignal=TRUE,
        	BGMedianSignal=TRUE,
        	BGUsed=TRUE,
		targets.micro)
graphics.off()

}
}
\keyword{documentation}
\keyword{utilities}

