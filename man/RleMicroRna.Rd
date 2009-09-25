\name{RleMicroRna}
\alias{RleMicroRna}
\title{ Relative Log Expression }
\description{
  RLE: Relative Log Expression 
}
\usage{
RleMicroRna(object, maintitle, colorfill)
}
\arguments{
  \item{object}{An expression matrix  }
  \item{maintitle}{title of the plot  }
  \item{colorfill}{color of the plot  }
}
\details{
	Each Boxplot corresponds to a sample and displays
	the Relative Log Expression computed for every spot
	in the array as the difference between the spot intensity and
	the median intensity for the same feature accros all the arrays. Since 
	majority of the spots are expected not to be differentially
	expressed, the plot shoud show boxplots centered around zero
	and all of them having the aproximately the same dispersion. 
	An array showing greater dispersion than the other, or 
	being not centered at zero could have quality problems. 
}
\references{  
	Boldstad B.M., Collin F., Brettschneider J., Simpson, K.,
	Cope L., Irizarry R. A., Speed T. P. Quality Assesement of
	Affymetrix GeneChip Data. In Bioinformatics and Computational
	Biology Solutions Using R and Bioconductor. (eds.)  Gentleman R., 
	Carey V. J., Huber W., Irizarry R. A., Dudoit S. (2005). Springer. 
}
\author{ Pedro Lopez-Romero }
\examples{
		data(dd.micro)
                RleMicroRna(log2(dd.micro$Rb),
			maintitle="log2 Mean Signal RLE",
			colorfill="orange")
}
\keyword{documentation}
\keyword{utilities}

