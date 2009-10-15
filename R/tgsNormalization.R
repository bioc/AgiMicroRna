`tgsNormalization` <-
function(ddTGS,NORMmethod="quantile",
		makePLOTpre=FALSE,makePLOTpost=FALSE,targets,verbose=FALSE){


	if (!is(ddTGS, "RGList")){
	  stop("'input' must be a RGList")
   	 	if (is.null(dim(ddTGS)[1])) {
        		stop("'input' is empty")
	 	}
	}
	if(NORMmethod != "none" && NORMmethod != "quantile" && NORMmethod != "scale"){
		stop("NORMmethod should be one of 'none', 'quantile','scale'" )	
	}
	
	if(!missing(makePLOTpre)) {
	if(makePLOTpre){
	
	MMM=log2(ddTGS$G)
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
		exprsNORM=normalizeBetweenArrays(ddTGS$G,method=NORMmethod) 
			if(NORMmethod == "quantile"){
				ddNORM$G=log2(exprsNORM)
			}else{
				ddNORM$G=exprsNORM
			}
		rm(exprsNORM)

	}else{
		ddNORM=ddTGS
		ddNORM$G=log2(ddTGS$G)
		 
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
	MMM=ddNORM$G

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

  ddNORM$Gb=ddNORM$G
  ddNORM$R=ddNORM$G 
  ddNORM$Rb=ddNORM$G
return(ddNORM)

} # end function 

