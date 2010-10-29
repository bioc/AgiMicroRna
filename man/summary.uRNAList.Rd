\name{summary.uRNAList}
\alias{summary.uRNAList}
\title{Summaries of Microarray Data Objects}
\description{
Briefly summarize microarray data objects.
}
\usage{
\method{summary}{uRNAList}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{uRNAList}}
  \item{...}{other arguments are not used}
}
\details{
The data objects are summarized as if they were lists, i.e., brief information about the length and type of the components is given.
This function and this file, has been borrowed from the files
created by Gordon Smyth for the limma package. 
}
\value{
A table.
}
\author{Pedro Lopez-Romero}
\keyword{methods}
