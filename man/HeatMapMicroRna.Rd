\name{HeatMapMicroRna}
\alias{HeatMapMicroRna}
\title{ HeatMap  }
\description{ Creates a HeatMap graph using the 'heatmap.2' function
}
\usage{
HeatMapMicroRna(object, size, maintitle)
}
\arguments{
  \item{object}{ A expression Matrix}
  \item{size}{ number of highest variance genes to be considered in the plot }
  \item{maintitle}{title of the plot  }
}
\author{ Pedro Lopez-Romero }
\seealso{ heatmap.2 }
\examples{
	data(ddPROC)
        HeatMapMicroRna(ddPROC$TGS,
		size=100,
        	maintitle="100 High Var genes")
}
\keyword{documentation}
\keyword{utilities}

