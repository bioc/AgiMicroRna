\name{mvaBASIC}
\alias{mvaBASIC}
\title{ MVA plot }
\description{
	For each array, the M value is computed for every spot as the difference
	between the spot intensity in the array and the averaged intensity for that feature over the whole
	set of arrays. It does not make a distintion between the different kind of features in the array 
	as the mvaMicroRna() does. 
}
\usage{
mvaBASIC(object, colorfill, maintitle)
}
\arguments{
  \item{object}{An expression matrix in log2 scale }
  \item{colorfill}{ color of the plot }
  \item{maintitle}{ title of the plot }
}
\author{ Pedro Lopez-Romero }
\examples{
	data(dd.micro)
	op=par(mfrow=c(1,1),ask=TRUE)
	mvaBASIC(log2(dd.micro$Rb),
	colorfill="red",
        maintitle=" log2 Mean Signal")
	par(op)
}
\keyword{documentation}
\keyword{utilities}
