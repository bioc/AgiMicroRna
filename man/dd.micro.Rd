\name{dd.micro}
\alias{dd.micro}
\docType{data}
\title{ data example (uRNAList) }
\description{
	Data, extracted from scanned images using Agilent Feature Extraction Software, are stored
	in a uRNAList object. 
}
\usage{data(dd.micro)}
\details{
	A data example is provided. The data example includes 3 experimental
	conditions with two replicates. 

  For these data, chips were scanned using the Agilent G2567AA
  Microarray Scanner System (Agilent Technologies)
  Image analysis and data collection were carried out
  using the Agilent Feature Extraction 9.1.3.1. (AFE).

	Data, colected with the Agilent Feature Extraction Software, are stored
	in a uRNAList object with the following components: 
\describe{
  \item{uRNAList\$TGS }{matrix, 'gTotalGeneSignal'}
  \item{uRNAList\$TPS }{matrix, 'gTotalProbeSignal'}
  \item{uRNAList\$meanS}{matrix, 'gMeanSignal'}
  \item{uRNAList\$procS}{matrix, 'gProcessedSignal'}
  \item{uRNAList\$targets }{data.frame, 'FileName'}
  \item{uRNAList\$genes\$ProbeName }{vector of characters, 'AGilent Probe Name'}
  \item{uRNAList\$genes\$GeneName }{vector of characters, 'microRNA Name'}
  \item{uRNAList\$genes\$ControlType }{vector of integers, '0'= Feature, '1'= Positive control, '-1'= Negative control}
  \item{uRNAList\$other\$gIsGeneDetected }{matrix, FLAG to classify signal if 'IsGeneDetected=1' or 'not=0'}
  \item{uRNAList\$other\$gIsSaturated }{matrix, FLAG to classify signal if 'IsSaturated = 1' or 'not=0'}
  \item{uRNAList\$other\$gIsFeatPopnOL }{matrix, FLAG to classify signal if 'IsFeatPopnOL = 0' or 'not=1'}
  \item{uRNAList\$other\$gIsFeatNonUnifOL }{matrix, FLAG to classify signal if 'gIsFeatNonUnifOL = 0' or 'not=1'}
  \item{uRNAList\$other\$gBGMedianSignal }{matrix, gBGMedianSignal }
  \item{uRNAList\$other\$gBGUsed }{matrix, gBGUsed }
}
}
\author{ Pedro Lopez-Romero }
\seealso{
	\code{readMicroRnaAFE.Rd}
}
\keyword{datasets}
