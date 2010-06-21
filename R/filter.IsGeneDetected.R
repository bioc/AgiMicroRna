`filter.IsGeneDetected` <-
function(ddFILT,limIsGeneDetected,targets,verbose,writeout){

	 
	minFLAGisf=1 	# gIsFound: FLAG ok: 1 = feature FOUND (58)

	if (!is(ddFILT, "RGList")){
	  stop("'input' must be a RGList")
   	 	if (is.null(dim(ddFILT)[1])) {
        		stop("'input' is empty")
	 	}
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

	FLAG=ddFILT$other$gIsGeneDetected
  	indexSNR=apply(FLAG,1,filterFLAG.micro,GErep,nGE,minFLAGisf,limIsGeneDetected)
  	selSNR=which(unlist(indexSNR) == 1) # posiciones que pasan filtro 

# WRITING OUT THE REMOVED DATA only geneIDs and selected flags

  if(writeout){
	 outfile="IsNOTGeneDetected.txt" 
	 write.filt.out.miRNA(ddFILT,selSNR,outfile,FLAG,targets) 
 }
	if(verbose){
	cat("FILTERING BY IsGeneDetected FLAG","\n") 
	cat("\n")
	cat("	FLAG FILTERING OPTIONS - FLAG OK = 1 - limIsGeneDetected: ",limIsGeneDetected,"%","\n")
	cat("	FEATURES AFTER IsGeneDetected FILTERING: ",sum(indexSNR),"\n")
	cat("	NON Gene Detected :",length(ddFILT$genes$ProbeName[-selSNR]),"\n")
	cat("------------------------------------------------------","\n")
	}

# EXTRACTING THE SELECTED PROBES  
	
ddFILT=ddFILT[selSNR,]	 
return(ddFILT)

}  # end function 

