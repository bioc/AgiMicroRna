\name{pvalHistogram}
\alias{pvalHistogram}
\title{ Histogram of the p values }
\description{
	Creates an histogram of the pvalues. For multiple contrats, creates
	an histogram for every t.test pvalue (separate) or a single histogram
	for the F.test pvalue (nestedF).
	A uniform histogram will indicate no differential expression in the
	data set, whereas a right skewed histogram, will indicate some significant
	differential expression
}
\usage{
pvalHistogram(fit2, DE, PVcut, DEmethod, MTestmethod, CM,verbose=FALSE)
}
\arguments{
  \item{fit2}{ MArrayLM object}
  \item{DE}{ TestResults object }
  \item{PVcut}{  limit p value to declare significant features}
  \item{DEmethod}{ method for decideTests, only 'separate' or 'nestedF' are implemented }
  \item{MTestmethod}{ method for multiple test }
  \item{CM}{ contrast matrix }
  \item{verbose}{logical, if \code{TRUE} prints out output} 
}
\author{ Pedro Lopez-Romero }
\seealso{An overview of miRNA differential expression analysis is given in
        \code{basicLimma}
	An example of how to get the 'TestResults' object is in
	\code{getDecideTests} 
}
\examples{
\dontrun{
pvalHistogram(fit2,DE,PVcut=0.10,
     DEmethod="separate",MTestmethod="BH",CM)
}
}
\keyword{documentation}
\keyword{utilities}

