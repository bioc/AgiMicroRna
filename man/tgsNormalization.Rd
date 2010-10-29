\name{tgsNormalization}
\alias{tgsNormalization}
\title{ Normalization Between Arrays }
\description{
  Normalization between arrays of the Total Gene Signal. The function
  is a wrapper of the 'limma' 'normalizeBetweenArrays' with ('none','quantile','scale')
  methods 
}
\usage{
tgsNormalization(ddTGS, NORMmethod = "quantile", makePLOTpre = FALSE, makePLOTpost = FALSE, targets,verbose=FALSE)
}
\arguments{
  \item{ddTGS}{uRNAList, containing the output from \code{tgsMicroRna} }
  \item{NORMmethod}{ character specifying the normalization method, 'none','quantile','scale'. The 
			default is \code{quantile} }
  \item{makePLOTpre}{ logical, if \code{TRUE} QC plots with the Raw Total Gene Signal are displayed  }
  \item{makePLOTpost}{ logical, if \code{TRUE} QC plots with the Normalized Total Gene Signal are displayed  }
  \item{targets}{ data.frame with the target structure }
  \item{verbose}{logical, if \code{TRUE} prints out output}
}

\value{
	A uRNAList object containing the Normalized Total Gene Signal in log 2 scale 
}
\references{ 
	Smyth, G. K. (2005). Limma: linear models for microarray data. In:
  	'Bioinformatics and Computational Biology Solutions Using R and
 	Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry,
  	W. Huber (eds), Springer, New York, pages 397 - 420

	Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data.
	Methods 31, 265-273.
}
\author{ Pedro Lopez-Romero }
\examples{
\dontrun{
data(dd.micro)
data(targets.micro)
ddTGS=tgsMicroRna(dd.micro,half=TRUE,makePLOT=FALSE,verbose=FALSE)

ddNORM=tgsNormalization(ddTGS,'quantile',
                makePLOTpre=FALSE,makePLOTpost=TRUE,targets.micro,verbose=TRUE)
graphics.off()

}
}
\keyword{documentation}
\keyword{utilities}

