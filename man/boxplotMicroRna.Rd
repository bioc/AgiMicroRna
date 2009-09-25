\name{boxplotMicroRna}
\alias{boxplotMicroRna}
\title{ Boxplot }
\description{
	It creates a Boxplot using the  matrix columns as input
}
\usage{
boxplotMicroRna(object, maintitle, colorfill, xlab, ylab)
}
\arguments{
  \item{object}{A matrix containing by columns the expression arrays in log2 scale}
  \item{maintitle}{character to indicate the title of the graph  }
  \item{colorfill}{ color to fill the boxplot }
  \item{xlab}{ title for the x axe }
  \item{ylab}{ title for the y axe }
}
\author{ Pedro Lopez-Romero }
\examples{
		data(dd.micro)
		MMM=log2(dd.micro$Rb)
                boxplotMicroRna(MMM,
                	maintitle="log2 Mean Signal",
                	colorfill="green",
			xlab="Samples",
			ylab="expression")
}
\keyword{documentation}
\keyword{utilities}
