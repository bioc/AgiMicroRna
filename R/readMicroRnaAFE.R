`readMicroRnaAFE` <-
function(targets,verbose=FALSE){

	if (!is(targets, "data.frame")){
        stop("'targets' must be a data.frame")
	}

dd=read.agiMicroRna(targets,
             columns=list(TGS="gTotalGeneSignal",
			         TPS="gTotalProbeSignal",
			         meanS="gMeanSignal",
			         procS="gProcessedSignal"),
	           other.columns=list(IsGeneDetected="gIsGeneDetected",
				        IsSaturated="gIsSaturated",
				        IsFeatNonUnifOF="gIsFeatNonUnifOL",
				        IsFeatPopnOL="gIsFeatPopnOL",
				        BGKmd="gBGMedianSignal",
				        BGKus="gBGUsed"),
	             annotation = c( "ControlType", "ProbeName","GeneName"),
	             verbose=TRUE)

  if(verbose){
		  cat("","\n")
		cat("  uRNAList:","\n")
		cat("	dd$TGS:		'gTotalGeneSignal' ","\n")
		cat("	dd$TPS:		'gTotalProbeSignal' ","\n")
		cat("	dd$meanS:	'gMeanSignal' ","\n")
		cat("	dd$procS:	'gProcessedSignal' ","\n")
		cat("","\n")
		}
return(dd)

} # end function 

