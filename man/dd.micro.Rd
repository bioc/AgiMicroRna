\name{dd.micro}
\alias{dd.micro}
\docType{data}
\title{ data example (RGList) }
\description{
	Data, extracted from scanned images using Agilent Feature Extraction Software, are stored
	in a RGList object. 
}
\usage{data(dd.micro)}
\details{
	A data example is provided. The data example includes 3 experimental
	conditions with two replicates. 

  Chips were scanned using the Agilent G2567AA
  Microarray Scanner System (Agilent Technologies)
  Image analysis and data collection were carried out
  using the Agilent Feature Extraction 9.1.3.1. (AFE).

	Data, colected with the Agilent Feature Extraction Software, are stored
	in a RGList object with the following components: 

   - dd.micro\$R: 		       		 'gTotalGeneSignal'  
   - dd.micro\$G:           		         'gTotalProbeSignal'  
   - dd.micro\$Rb:          		         'gMeanSignal'  
   - dd.micro\$Gb:          		         'gProcessedSignal' 
   - dd.micro\$targets       		         'targets'
   - dd.micro\$genes\$ProbeName		     	 'Probe Name'
   - dd.micro\$genes\$GeneName		         'microRNA Name'
   - dd.micro\$genes\$ControlType		 'FLAG to specify the sort of feature'
   - dd.micro\$other\$gIsGeneDetected	         'FLAG IsGeneDetected'
   - dd.micro\$other\$gIsSaturated		 'FLAG IsSaturated'
   - dd.micro\$other\$gIsFeatNonUnifOL	         'FLAG IsFeatNonUnifOL'
   - dd.micro\$other\$gIsFeatPopnOL	         'FLAG IsFeatPopnOL'
   - dd.micro\$other\$chr\_coord		 'CHR coordinate (obtained from Agilent data files)
   - dd.micro\$other\$gBGMedianSignal	         'gBGMedianSignal'
   - dd.micro\$other\$gBGUsed	                 'gBGUsed'
}
\author{ Pedro Lopez-Romero }
\seealso{
	\code{readMicroRnaAFE.Rd}
}
\keyword{datasets}
