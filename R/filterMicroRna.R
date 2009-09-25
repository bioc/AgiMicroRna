`filterMicroRna` <-
function(ddNORM,
	dd,
	control=TRUE,
	IsGeneDetected=TRUE,
	wellaboveNEG=FALSE,
	limIsGeneDetected=75,
	limNEG=75,
	makePLOT=FALSE,
	targets,
	verbose=FALSE) {

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

	if (!is(ddNORM, "RGList")){
	  stop("'input' must be a RGList")
   	 	if (is.null(dim(ddNORM)[1])) {
        		stop("'input' is empty")
	 	}
	}


# FILTERING STEPS 

	if(!missing(control)){
	if(control){
		ddFILT=filter.control.miRNA(ddNORM,targets,verbose)	
			} # FILTERct 
			}

	if(!missing(IsGeneDetected)){
	if(IsGeneDetected){
		ddFILT=filter.IsGeneDetected(ddFILT,limIsGeneDetected,targets,verbose)
			} # wellaboveBG 
			}

	if(!missing(wellaboveNEG)){
	if(wellaboveNEG){
		ddFILT=filter.wellaboveNEG.miRNA(ddFILT,dd,limNEG,SDtimes=1.5,targets,verbose)
			} # wellaboveNEG
			}

if(!missing(makePLOT)) {
	if(makePLOT){
	

	colorfill="yellow"		
	maintitle="FILTERED SIGNAL"
	
		dev.new()
		plotDensityMicroRna(ddFILT$G,maintitle)

		dev.new()
		boxplotMicroRna(ddFILT$G,maintitle,colorfill)
	
		dev.new()
		hierclusMicroRna(ddFILT$G,targets$GErep,methdis="euclidean",
        	methclu="complete",sel=FALSE,100)

		maintitle="FILTERED SIGNAL - RLE "
		dev.new()
		RleMicroRna(ddFILT$G,maintitle,colorfill)
	}
	}
	return(ddFILT) 

} # end function 

