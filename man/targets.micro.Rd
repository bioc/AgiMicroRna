\name{targets.micro}
\alias{targets.micro}
\docType{data}
\title{ Example of target file}
\description{ Example of target file 
}
\usage{data(targets.micro)}
\format{
  A data frame with 4 observations on the following 5 variables.
  \describe{
    \item{\code{FileName}}{names of the Files \code{Ast.txt} \code{Bst.txt} \code{Aunst.txt} \code{Bunst.txt}}
    \item{\code{Treatment}}{Assigns level for Treatment Effect to each File (mandatory)}
    \item{\code{GErep}}{a numeric vector tha numerates the FACTOR of the Treatment Effect (mandatory)}
    \item{\code{Subject}}{Assigns level for Subject Effect to each File}
  }
}
\details{
	It is a tab-delimited text format file. 
	The target file is created by the user with the intention of
	carrying out a differential expression analysis  
	in future steps using 'limma'. Here is where the 
	factors that are going to be included in the 
	linear model that is fitted to each gen are specified. 
	The targets file assigns each  data file to  a particular 
	experimental conditions. 
	First column 'FileName' is mandatory and includes the
	image data files names. Second column
	'Treatment' is also mandatory and includes the treament
	effect. Third column 'GErep' is also mandatory, and
	includes the treatment effect in a numeric code, from 
	1 to n, being n the number of Treatment effect levels.
}
\references{
	Gordon K. Smyth, M. Ritchie, N. Thorne, J. Wettenhall (2007). limma:
	Linear Models for Microarray Data User's Guide. 
}
\author{ Pedro Lopez-Romero }
\seealso{
        \code{readTargets}
}

\keyword{datasets}
