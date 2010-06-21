`filter.wellaboveNEG.miRNA` <-
function(ddFILT,dd,limNEG,SDtimes,targets,verbose,writeout){

	 
	if (!is(ddFILT, "RGList")){
	  stop("'input ddFILT' must be a RGList")
   	 	if (is.null(dim(ddFILT)[1])) {
        		stop("'input' is empty")
	 	}
	}

	if (!is(dd, "RGList")){
	  stop("'input dd' must be a RGList")
   	 	if (is.null(dim(dd)[1])) {
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

	indexneg=which(dd$genes$ControlType == -1)
	NEG=log2(dd$Rb[indexneg,])  # dd$G: MeanSignal 
	MeanNeg=apply(NEG,2,mean)
	SdNeg=apply(NEG,2,sd)
	Limit=MeanNeg + SDtimes*(SdNeg)
  
  	
  	FLAG=ddFILT$G 
  	indexSNR=apply(FLAG,1,filterWellAboveSIGNALv2,GErep,nGE,Limit,limNEG)
  	selSNR=which(unlist(indexSNR) == 1) #posiciones que pasan filtro 

		
# WRITING OUT THE WellAboveNeg.out filtered DATA 

if(writeout){
	outfile="IsNOTWellAboveNEG.txt"
	VALUES=round(ddFILT$G,3)
	write.filt.out.miRNA(ddFILT,selSNR,outfile,VALUES,targets) 
}

	      
		if(verbose){
		cat("FILTERING BY WellAboveNeg filterWellAboveSIGNALv2 ~ FLAG","\n")
		cat("\n")
		cat("	FLAG FILTERING OPTIONS - limNEG: ",limNEG,"%","\n")
  		cat("	Limit computed as MeanNeg + ",SDtimes," x (SDNeg) ","\n")
  		
  		cat("	Limit in each array: ",round(Limit,2),"\n")
  		cat("\n")
		cat("PROBES AFTER WellAboveNeg FILTERING: ",sum(indexSNR),"\n")
		cat("WellAboveNeg OUT :",length(ddFILT$genes$ProbeName[-selSNR]),"\n")
		cat("------------------------------------------------------","\n")
		}

# EXTRACTING THE SELECTED PROBES  
	
	ddFILT=ddFILT[selSNR,]
	return(ddFILT)

} # end function 

