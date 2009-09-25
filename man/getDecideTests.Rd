\name{getDecideTests}
\alias{getDecideTests}
\title{Differential expression analysis an multiplicity of the tests }
\description{
	It Uses the \code{decideTests} function of the 'limma' package to
	classify the list of genes as up, down or not significant after 
	correcting by the multiplicity of the tests. 
}
\usage{
getDecideTests(fit2, DEmethod, MTestmethod, PVcut,verbose=FALSE)
}
\arguments{
  \item{fit2}{MArrayLM object }
  \item{DEmethod}{method for \code{decideTests}, only 'separate' or 'nestedF' are implemented.
	see \code{decideTests} in limma package.  
	} 
  \item{MTestmethod}{method for multiple test, choices  are 'none','BH', 'BY', ...
	see \code{p.adjust}
	}
  \item{PVcut}{p value threshold to declare significant features }
  \item{verbose}{logical, if \code{TRUE} prints out output}  
  
}
\value{
	A 'TestResults' object of the 'limma' package 
	It prints out the number of UP and DOWN genes for every contrasts
	according to the p value limit specified
}
\references{  
	Smyth, G. K. (2005). Limma: linear models for microarray data. In:
        'Bioinformatics and Computational Biology Solutions using R and
        Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.
        Huber (eds), Springer, New York, pages 397--420.
}
\author{ Pedro Lopez-Romero }
\seealso{An overview of miRNA differential expression analysis is given in
	\code{basicLimma}
}
\examples{
\dontrun{
DE=getDecideTests(fit2,
        DEmethod="separate",
        MTestmethod="BH",
        PVcut=0.10,
	verbose=TRUE)
}
}
\keyword{documentation}
\keyword{utilities}

