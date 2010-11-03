\name{uRNAList}
\docType{class}
\alias{uRNAList}
\alias{uRNAList-class}
\alias{[,uRNAList-method}
\alias{show,uRNAList-method}

\title{uRNAList - class}

\description{
A list-based class (similar to the \code{RGList} class in \code{limma} package) for storing microRNA data
of Agilent chips 
\code{uRNAList} objects are created by \code{read.agiMicroRna}
}

\section{uRNAList Components}{
\code{uRNAList} objects are created by \code{new("uRNAList",Newagi)} where \code{Newagi} is a list. 
  with the following components:
\describe{
  \item{uRNAList\$TGS}{matrix, 'gTotalGeneSignal'}
  \item{uRNAList\$TPS}{matrix, 'gTotalProbeSignal'}
  \item{uRNAList\$meanS}{matrix, 'gMeanSignal'}
  \item{uRNAList\$procS}{matrix, 'gProcessedSignal'}
  \item{uRNAList\$targets}{data.frame, 'FileName'}
  \item{uRNAList\$genes\$ProbeName}{vector of characters, 'AGilent Probe Name'}
  \item{uRNAList\$genes\$GeneName}{vector of characters, 'microRNA Name'}
  \item{uRNAList\$genes\$ControlType}{vector of integers, '0'= Feature, '1'= Positive control, '-1'= Negative control}
  \item{uRNAList\$other\$gIsGeneDetected}{matrix, FLAG to classify signal if 'IsGeneDetected=1' or 'not=0'}
  \item{uRNAList\$other\$gIsSaturated}{matrix, FLAG to classify signal if 'IsSaturated = 1' or 'not=0'}
  \item{uRNAList\$other\$gIsFeatPopnOL}{matrix, FLAG to classify signal if 'IsFeatPopnOL = 0' or 'not=1'}
  \item{uRNAList\$other\$gIsFeatNonUnifOL}{matrix, FLAG to classify signal if 'gIsFeatNonUnifOL = 0' or 'not=1'}
  \item{uRNAList\$other\$gBGMedianSignal}{matrix, gBGMedianSignal }
  \item{uRNAList\$other\$gBGUsed}{matrix, gBGUsed }
}

}

\author{Pedro Lopez-Romero}
\examples{
\dontrun{
    data(dd.micro)
}
}
\keyword{classes}
