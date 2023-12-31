`filter.control.miRNA` <-
function(ddNORM,targets,verbose,writeout){

	if (!is(ddNORM, "uRNAList")){
	  stop("'input' must be a uRNAList")
   	 	if (is.null(dim(ddNORM)[1])) {
        		stop("'input' is empty")
	 	}
	}
	
	if(missing(targets)){
		stop("'targets' is missing ")
	}

	if("GErep" %in% colnames(targets)){
		g1=targets$GErep  # g1 must be numeric, from 1:n
		g2=rownames(targets)
	}else{
		stop("'targets' needs 'GErep' field")
	}

	cat("\n")
	cat("FILTERING BY ControlType","\n")
  
	selSNR=which(ddNORM$genes$ControlType==0)
 
	ddFILT=ddNORM[selSNR,]

		
# WRITING OUT THE RAW DATA WITHOUT CONTROLS 

if(writeout){
	write.control.out.miRNA(ddFILT,selSNR,targets)
}	
		if(verbose){
    		cat("\n")
		cat("   FEATURES BEFORE FILTERING: ",dim(ddNORM)[1],"\n")
		cat(" 	FEATURES AFTER ControlType FILTERING: ",dim(ddFILT)[1],"\n")
		cat("------------------------------------------------------","\n")
		}

	return(ddFILT)
} # end function 

