`tgsMicroRna` <-
function(dd,offset=0,half=TRUE,makePLOT=FALSE,verbose=FALSE){

	if (!is(dd, "RGList")){
	  stop("'input' must be a RGList")
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
# cat("		both ddTGS$G & ddTGS$R contains gTotalGeneSignal","\n")
 }
 
ug=which(duplicated(dd$genes$GeneName)==FALSE)
nGEN=length(ug)
ddTGS=dd[ug,] # RGlist with only TotalGeneSignal 

if(half){
	# cat("ddTGS signal with 'half method'","\n")

	for(i in 1:dim(ddTGS)[2]){
	index=which(ddTGS$R[,i] < 0.5)
	ddTGS$R[index,i]=0.5 
	} 
}else{
	# cat("ddTGS signal with offset of 'abs(min(ddTGS$R))' + ",offset,"\n")
	min=min(ddTGS$R)
	for(i in 1:dim(ddTGS)[2]){
	ddTGS$R[,i]=ddTGS$R[,i]+(abs(min)+offset)
	} 
}

ddTGS$G=ddTGS$R 
ddTGS$R=ddTGS$R
ddTGS$Gb=ddTGS$R 
ddTGS$Rb=ddTGS$R

nARR=dim(ddTGS)[2]
geneNames=list(c(dd$genes$GeneName[ug]),c(1:nARR)) 	# GENES: condensated in one final gene interrogated by multiple probes 	

if(!missing(makePLOT)) {
	  if(makePLOT){
		
	MMM=log2(ddTGS$R)
	maintitle="TotalGeneSignal"
	colorfill="green"

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		dev.new()
		ddaux=ddTGS
		ddaux$G=MMM
		mvaMicroRna(ddaux,maintitle,TRUE)
		rm(ddaux)

	}
}

return(ddTGS)

} # end function 

