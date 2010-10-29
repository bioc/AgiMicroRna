`tgsNormalization` <-
function(ddTGS,NORMmethod="quantile",
		makePLOTpre=FALSE,makePLOTpost=FALSE,targets,verbose=FALSE){


	if (!is(ddTGS, "uRNAList")){
	  stop("'input' must be a uRNAList")
   	 	if (is.null(dim(ddTGS)[1])) {
        		stop("'input' is empty")
	 	}
	}
	if(NORMmethod != "none" && NORMmethod != "quantile" && NORMmethod != "scale"){
		stop("NORMmethod should be one of 'none', 'quantile','scale'" )	
	}
	
	if(!missing(makePLOTpre)) {
	if(makePLOTpre){
	
	MMM=log2(ddTGS$TGS)
	colorfill="blue"

	maintitle="NOT NORM."
		
		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		mvaBASIC(MMM,colorfill,maintitle) 

		dev.new()
		hierclusMicroRna(MMM,targets$GErep,methdis="euclidean",
        	methclu="complete",sel=FALSE,100)

		maintitle="NOT NORM - RLE "
		dev.new()
		RleMicroRna(MMM,maintitle,colorfill)

		rm(MMM)
	}
	}

	if(NORMmethod != "none"){
		ddNORM=ddTGS
		exprsNORM=normalizeBetweenArrays(ddTGS$TGS,method=NORMmethod) 
			if(NORMmethod == "quantile"){
				ddNORM$TGS=log2(exprsNORM)
			}else{
				ddNORM$TGS=exprsNORM
			}
		rm(exprsNORM)

	}else{
		ddNORM=ddTGS
		ddNORM$TGS=log2(ddTGS$TGS)
		 
	}

	if(verbose){
	cat("------------------------------------------------------","\n")
	cat("	NORMMALIZATION:	",NORMmethod,"\n")
	cat("	OUTPUT in log-2 scale","\n")
	cat("------------------------------------------------------","\n")
        }

	if(!missing(makePLOTpost)) {
	if(makePLOTpost){
	colorfill="red"

	maintitle="NORMALIZED SIGNAL"
	MMM=ddNORM$TGS

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		mvaBASIC(MMM,colorfill,maintitle) 

		dev.new()
		hierclusMicroRna(MMM,targets$GErep,methdis="euclidean",
        	methclu="complete",sel=FALSE,100)

		maintitle="NORM DATA - RLE "
		dev.new()
		RleMicroRna(MMM,maintitle,colorfill)

		rm(MMM)
		
	}
	}

  ddNORM$TGS=ddNORM$TGS
  ddNORM$TGS=ddNORM$TGS 
  ddNORM$TGS=ddNORM$TGS 
return(ddNORM)

} # end function 

