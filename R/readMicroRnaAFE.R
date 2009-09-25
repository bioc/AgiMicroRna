`readMicroRnaAFE` <-
function(targets,verbose=FALSE){

	if (!is(targets, "data.frame")){
        stop("'targets' must be a data.frame")
	}

  ddaux=read.maimages(files=targets$FileName,source="agilent",
	           other.columns=list(IsGeneDetected="gIsGeneDetected",
				        IsSaturated="gIsSaturated",
				        IsFeatNonUnifOF="gIsFeatNonUnifOL",
				        IsFeatPopnOL="gIsFeatPopnOL",
				        ChrCoord="chr_coord",
				        BGKmd="gBGMedianSignal",
				        BGKus="gBGUsed"),
	         columns=list(Rf="gTotalGeneSignal",
			         Gf="gTotalProbeSignal",
			         Rb="gMeanSignal",
			         Gb="gProcessedSignal"),
	         verbose=TRUE,sep="\t",quote="")

        dd=new("RGList")

                dd$R=ddaux$R
                dd$G=ddaux$G
                dd$Rb=ddaux$Rb
                dd$Gb=ddaux$Gb

                dd$targets=ddaux$targets
                dd$genes=ddaux$genes[,c(4,5,6)]
                dd$other=ddaux$other

        rm(ddaux)

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

