`tgsMicroRna` <-
function(dd,offset=0,half=TRUE,makePLOT=FALSE,verbose=FALSE){

	if (!is(dd, "uRNAList")){
	  stop("'input' must be a uRNAList")
   	 	if (is.null(dim(dd)[1])) {
        		stop("'input' is empty")
	 	}
	}

TT=dim(dd)[1]
uniqueProbe=unique(dd$genes$ProbeName)
LUP=length(uniqueProbe)
uniqueGene=unique(dd$genes$GeneName)
LUG=length(uniqueGene)

if(verbose){
 cat("\n")
 cat("GETTING Agilent Feature Extraction TotalGeneSignal","\n")
 cat("\n")
 cat("	Total Probes:	",TT,"\n")
 cat("	Unique Probe: 	",LUP,"\n")
 cat("	Unique Gene: 	",LUG,"\n")
 cat("\n")
 }
 
ug=which(duplicated(dd$genes$GeneName)==FALSE)
nGEN=length(ug)
ddTGS=dd[ug,] # uRNAList with only TotalGeneSignal 

if(half){
	# cat("ddTGS signal with 'half method'","\n")

	for(i in 1:dim(ddTGS)[2]){
	index=which(ddTGS$TGS[,i] < 0.5)
	ddTGS$TGS[index,i]=0.5 
	} 
}else{
	# cat("ddTGS signal with offset of 'abs(min(ddTGS$TGS))' + ",offset,"\n")
	min=min(ddTGS$TGS)
	for(i in 1:dim(ddTGS)[2]){
	ddTGS$TGS[,i]=ddTGS$TGS[,i]+(abs(min)+offset)
	} 
}

ddTGS$TGS= ddTGS$TGS 
ddTGS$TGS= ddTGS$TGS
ddTGS$TGS= ddTGS$TGS
ddTGS$TGS= ddTGS$TGS

nARR=dim(ddTGS)[2]
geneNames=list(c(dd$genes$GeneName[ug]),c(1:nARR)) 	# GENES: condensated in one final gene interrogated by multiple probes 	

if(!missing(makePLOT)) {
	  if(makePLOT){
		
	MMM=log2(ddTGS$TGS)
	maintitle="TotalGeneSignal"
	colorfill="green"

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		dev.new()
		ddaux=ddTGS
		ddaux$meanS=MMM
		mvaMicroRna(ddaux,maintitle,TRUE)
		rm(ddaux)

	}
}

return(ddTGS)

} # end function 

