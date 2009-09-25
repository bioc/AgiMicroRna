`cvArray` <-
function(ddDUP,foreground=c("MeanSignal","ProcessedSignal"),targets,verbose=FALSE){


	if (!is(ddDUP, "RGList")){
	  stop("'input' must be a RGList",call. = FALSE)
   	 	if (is.null(dim(ddDUP)[1])) {
        		stop("'input' is empty",call. = FALSE)
	 	}
	}
	foreground <- match.arg(foreground)

	if(foreground != "ProcessedSignal" && foreground !="MeanSignal" ){
		if(verbose){
		cat("Foreground: ",foreground,"\n")
		cat("\n")
		}
		stop("'Foreground' must be either 'ProcessedSignal' or 'MeanSignal'",call.=FALSE)	
	}else{ 
		if(foreground == "ProcessedSignal"){
	  		ddDUP$G=ddDUP$Gb
			if(verbose){
	  		cat("Foreground: ProcessedSignal","\n")
			cat("\n")
			}
		}else{
			ddDUP$G=ddDUP$Rb
			if(verbose){
			cat("Foreground: MeanSignal","\n")
			cat("\n")
			}
		}
	}

		# ControlType: 0 = signal, +1 = C.pos, -1 = C.neg (56)

		if(verbose){
		cat("	FILTERING BY ControlType FLAG","\n")
		cat("\n")
		cat(" RAW DATA: 			",dim(ddDUP)[1],"\n")
		}

		selSNR=which(ddDUP$genes$ControlType==0)
 		ddDUP=ddDUP[selSNR,]
 		nGEN=dim(ddDUP)[1]
 		
		if(verbose){
		cat(" PROBES without CONTROLS: 	",nGEN,"\n")
		}
		extra=which(duplicated(ddDUP$genes$ProbeName)==TRUE)
		if(length(extra) == 0){
			stop("NOT DUPLICATED ProbeName in chip")
			} 

		uniqueProbe=unique(ddDUP$genes$ProbeName)
		LUP=length(uniqueProbe)
		uniqueGene=unique(ddDUP$genes$GeneName)
		LUG=length(uniqueGene)
		if(verbose){
		cat("----------------------------------","\n")
			cat("  (Non-CTRL) Unique Probe: ",LUP,"\n")
			cat("  (Non-CTRL) Unique Genes: ",LUG,"\n")
		cat("----------------------------------","\n")
		}

		reps=table(ddDUP$genes$ProbeName)
		# reps[]=reps[]+1
		t=table(reps)

		rN=names(reps)
		Lreps=length(rN)

		if(verbose){
		cat("DISTRIBUTION OF REPLICATED NonControl Probes","\n")
		print(t)
		cat("------------------------------------------------------","\n")
		}

		cv.array(ddDUP,rN,targets,"ProbeName",verbose)

		reps=table(ddDUP$genes$GeneName)
		# reps[]=reps[]+1
		t=table(reps)

		rN=names(reps)
		Lreps=length(rN)

		if(verbose){
		cat("DISTRIBUTION OF REPLICATED Noncontrol Genes","\n")
		print(t)
		cat("------------------------------------------------------","\n")
		}
		
#		cv.array(ddDUP,rN,targets,"GeneName",verbose)

} # end function 

