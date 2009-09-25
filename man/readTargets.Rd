\name{readTargets}
\alias{readTargets}
\title{read the target file}
\description{
	The target file is a txt file created by the user 
	where every input file (array, sample) is
	attached to a experimental condition
}
\usage{
	readTargets(infile,verbose=FALSE)
}
\arguments{
  \item{infile}{ name of the target file, for instance 'targets.micro.txt' }
   \item{verbose}{logical, if \code{TRUE} prints out output}
}
\details{
  In the 'target' file (see Table 1 in vignette) we specify the 
  experimental conditions under which the data have been generated. 
  The target file MUST contain the following mandatory columns:
		-FileName	: Name of the array data file
		-Treatment	: Treatment effect  
		-GErep		: Treatment effect  in numeric code, from 
      				  '1' to 'n', being 'n' the number of the levels
				   of the treatment effect
      
  Other explanatory variables specifying the experimental conditions might be
  also included. 
}
\value{
	A 'data.frame' containing by the columns specified in the input file
	targets.txt. This 'targets.txt' file must be created by the user.  
}
\references{  
        Smyth, G. K. (2005). Limma: linear models for microarray data. In:
        'Bioinformatics and Computational Biology Solutions using R and
        Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.
        Huber (eds), Springer, New York, pages 397--420.
}
\author{ Pedro Lopez-Romero }
\seealso{
       An example of a target file can be found in \code{targets.micro}
}

\keyword{documentation}
\keyword{utilities}

