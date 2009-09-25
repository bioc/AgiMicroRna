\name{plotDensityMicroRna}
\alias{plotDensityMicroRna}
\title{ Density Plots of Intensity Signals }
\description{
	Creates a density plot with the arrays intensities
}
\usage{
plotDensityMicroRna(object, maintitle)
}
\arguments{
  \item{object}{An expression matrix, in log2 scale  }
  \item{maintitle}{title of the plot }
}
\author{ Pedro Lopez-Romero }
\examples{
	data(dd.micro)
        plotDensity(log2(dd.micro$Rb),maintitle="log2 Mean Signal")
}
\keyword{documentation}
\keyword{utilities}

