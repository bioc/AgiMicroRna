\name{dim.uRNAList}
\alias{dim.uRNAList}
\alias{length.uRNAList}
\title{Retrieve the Dimensions of an uRNAList Object}
\description{
Retrieve the number of rows (genes) and columns (arrays) for an uRNAList object.
}
\usage{
\method{dim}{uRNAList}(x)
\method{length}{uRNAList}(x)
}
\arguments{
  \item{x}{an object of class \code{uRNAList}}
}
\details{
This function and this file, has been borrowed from the files
created by Gordon Smyth for the limma package. 
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (arrays).
}
\author{Pedro Lopez-Romero}
\keyword{array}
