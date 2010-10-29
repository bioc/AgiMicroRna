\name{dimnames.uRNAList}
\alias{dimnames.uRNAList}
\alias{dimnames<-.uRNAList}
\title{Retrieve the Dimension Names of an uRNAList Object}
\description{
Retrieve the dimension names of a microarray data object.
}
\usage{
\method{dimnames}{uRNAList}(x)
\method{dimnames}{uRNAList}(x) <- value
}
\arguments{
  \item{x}{an object of class \code{uRNAList}}
  \item{value}{a possible value for \code{dimnames(x)}}  
}
\details{
The dimension names of a microarray object are the same as those of the most important matrix component of that object.
A consequence is that \code{rownames} and \code{colnames} will work as expected.
This function and this file, has been borrowed from the files
created by Gordon Smyth for the limma package. 
}
\value{
Either \code{NULL} or a list of length 2.
}
\author{Pedro Lopez-Romero}
\keyword{array}
