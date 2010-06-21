`readMicroRnaAFE` <-
function(targets,verbose=FALSE){

	if (!is(targets, "data.frame")){
        stop("'targets' must be a data.frame")
	}

dd=read.maimages(files=targets$FileName,source="agilent",
             columns=list(Rf="gTotalGeneSignal",
			         Gf="gTotalProbeSignal",
			         Rb="gMeanSignal",
			         Gb="gProcessedSignal"),
	           other.columns=list(IsGeneDetected="gIsGeneDetected",
				        IsSaturated="gIsSaturated",
				        IsFeatNonUnifOF="gIsFeatNonUnifOL",
				        IsFeatPopnOL="gIsFeatPopnOL",
				        BGKmd="gBGMedianSignal",
				        BGKus="gBGUsed"),
	             annotation = c( "ControlType", "ProbeName","GeneName"),
	             verbose=TRUE,sep="\t",quote="")

  if(verbose){
		  cat("","\n")
		cat("  RGList:","\n")
		cat("	dd$R:		'gTotalGeneSignal' ","\n")
		cat("	dd$G:		'gTotalProbeSignal' ","\n")
		cat("	dd$Rb:		'gMeanSignal' ","\n")
		cat("	dd$Gb:		'gProcessedSignal' ","\n")
		cat("","\n")
		}
return(dd)

} # end function 

