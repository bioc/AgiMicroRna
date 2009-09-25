\name{hierclusMicroRna}
\alias{hierclusMicroRna}
\title{ Hierarchical clustering  }
\description{  Hierarchical cluster of samples using the 'hclust' function
}
\usage{
hierclusMicroRna(object, GErep, methdis, methclu,sel, size)
}
\arguments{
  \item{object}{ An expression Matrix }
  \item{GErep}{ Numerical vector that relates each sample with its experimental condition }
  \item{methdis}{the distance measure to be used. Options are 'euclidean' and 'pearson'.
		 see 'dist' function  }
  \item{methclu}{the agglomeration method to be used by the 'hclust' function  }
  \item{sel}{logical, if \code{TRUE} selects the 'size' highest variance genes for the plot  }
  \item{size}{selects the 'size' highest variance genes for the plot if 'sel=TRUE'  }
}
\author{ Pedro Lopez-Romero }
\seealso{ hclust,dist }
\examples{
 	data(targets.micro)
        data(ddPROC)

	hierclusMicroRna(ddPROC$Rb,GErep,
		methdis="euclidean",
        	methclu="complete",
		sel=FALSE,100)
}
\keyword{documentation}
\keyword{utilities}


