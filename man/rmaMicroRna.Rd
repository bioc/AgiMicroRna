\name{rmaMicroRna}
\alias{rmaMicroRna}
\title{ Getting the Total Gene Signal by RMA algorithm }
\description{
	The function creates an RGList containing the TotalGeneSignal computed
	by the RMA algorithm. This signal can be used for the statistical analysis.  
}
\usage{
rmaMicroRna(dd, normalize, background)
}
\arguments{
  \item{dd}{RGList, containing the output from \code{readMicroRnaAFE} }
  \item{normalize}{logical, if \code{TRUE} the signal is normalized between arrays
	using the 'quantile' method}
  \item{background}{logical, if \code{TRUE} the signal is background corrected 
	by fitting a normal + exponential convolution model to a vector of 
	observed intensities}
}
\details{
	The function creates an RGList output that contains in the RGList$G, RGList$R,
        RGList$Gb & RGList$Rb slots the Total Gene Signal (TGS) computed by the RMA algorithm. 
	The function uses the robust multiarray average (RMA) method from the 'affy' package. 
	RMA obtains an estimate of the expression measure for each gene using all the replicated probes
	for that gene. First, RMA obtains a background corrected intensity by fitting a normal + exponential
	convolution model to a vector of observed intensities. The normal part represents
	the background and the exponential part represents the signal intensities. Then the
	arrays are normalized using 'quantile' normalization. Finally,
	for each probe set that interrogates the same microRNA, RMA fits a linear model to the background-corrected,
	normalized and log2 transformed probe intensities. This model produces an estimate
	of the gene signal taking into account the probe effect. The model parameters estimates are obtained
	by median polish. The estimates of the gene expression signals are referred as RMA estimates.
	Normally, each microRNA is interrogated by 16 probes either using 2 different probes, each
	of them replicated 8 times, or using 4 differnt probes  replicated 4 times. 
	First, function 'rmaMicroRna' obtains a background corrected signal using the 'rma.background.correct' function
	of the package 'preprocessCore' , then the signal is normalized bewtween
	arrays using the 'limma' function 'normalizeBetweenArrays' with the 'quantile' method. 
	Then, the median of the replicated probes is obtained, leading to either 2 or 4 different measures for each gene.
	These measures correspond to different probes for the same genes that are summarized into a single
	RMA linear model described above.
}
\value{
	RGList containing the Total Gene Signal computed
        by the RMA algorithm in log 2 scale.  
}
\references{ 
	Irizarry, R., Hobbs,B., Collin,F., Beazer-Barclay,Y., Antonellis,K., 
	Scherf,U., Speed,T. (2003) Exploration, normalization, and summaries 
	of high density oligonucleotide array probe level data. Biostatistics. 
	4, 249-264

	Gautier, L., Cope, L., Bolstad, B. M., and Irizarry, R. A.(2004). affy---analysis 
	of Affymetrix GeneChip data at the probe level. Bioinformatics 20, 3, 307-315.

	Bolstad B. M. (). preprocessCore: A collection of pre-processing functions.
	R package version 1.4.0

	Smyth, G. K. (2005). Limma: linear models for microarray data. In:  'Bioinformatics 
	and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, 
	V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397 - 420
}
\author{ Pedro Lopez-Romero }
\examples{
data(dd.micro)
ddTGS.rma=rmaMicroRna(dd.micro, normalize=TRUE, background=TRUE)
dim(ddTGS.rma)
RleMicroRna(ddTGS.rma$G,"RLE TGS.rma","blue")
}

\keyword{documentation}
\keyword{utilities}
