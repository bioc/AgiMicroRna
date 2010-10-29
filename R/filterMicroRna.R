`filterMicroRna` <-
function(ddNORM,
	dd,
	control,
	IsGeneDetected,
	wellaboveNEG,
	limIsGeneDetected,
	limNEG,
	makePLOT,
	targets,
	verbose,
  	writeout) {

	if(verbose){
	cat("FILTERING PROBES BY FLAGS","\n")
	cat("\n")
	}
	if(missing(targets)){
		stop("'targets' is missing ")
	}
	if("GErep" %in% colnames(targets)){
		g1=targets$GErep  # g1 must be numeric, from 1:n
		g2=rownames(targets)
		GErep=targets$GErep 
		nGE=sum(table(table(GErep)))
	}else{
		stop("'targets' needs 'GErep' field")
	}

	if (!is(ddNORM, "uRNAList")){
	  stop("'input' must be a uRNAList")
   	 	if (is.null(dim(ddNORM)[1])) {
        		stop("'input' is empty")
	 	}
	}


# FILTERING STEPS 

	if(!missing(control)){
	if(control){
		ddFILT=filter.control.miRNA(ddNORM,targets,verbose,writeout)	
			} # FILTERct 
			}

	if(!missing(IsGeneDetected)){
	if(IsGeneDetected){
		ddFILT=filter.IsGeneDetected(ddFILT,limIsGeneDetected,targets,verbose,writeout)
			} # wellaboveBG 
			}

	if(!missing(wellaboveNEG)){
	if(wellaboveNEG){
		ddFILT=filter.wellaboveNEG.miRNA(ddFILT,dd,limNEG,SDtimes=1.5,targets,verbose,writeout)
			} # wellaboveNEG
			}

if(!missing(makePLOT)) {
	if(makePLOT){
	

	colorfill="yellow"		
	maintitle="FILTERED SIGNAL"
	
		dev.new()
		plotDensityMicroRna(ddFILT$TGS,maintitle)

		dev.new()
		boxplotMicroRna(ddFILT$TGS,maintitle,colorfill)
	
		dev.new()
		hierclusMicroRna(ddFILT$TGS,targets$GErep,methdis="euclidean",
        	methclu="complete",sel=FALSE,100)

		maintitle="FILTERED SIGNAL - RLE "
		dev.new()
		RleMicroRna(ddFILT$TGS,maintitle,colorfill)
	}
	}
	return(ddFILT) 

} # end function 

