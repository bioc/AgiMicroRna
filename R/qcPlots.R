`qcPlots` <-
function(dd,offset=5,MeanSignal=TRUE,ProcessedSignal=FALSE,
		TotalProbeSignal=FALSE,TotalGeneSignal=FALSE,
		BGMedianSignal=FALSE,BGUsed=FALSE,targets){

	if (!is(dd, "uRNAList")){
	  stop("'input' must be a uRNAList")
   	 	if (is.null(dim(dd)[1])) {
        		stop("'input' is empty")
	 	}
	}
	

if(MeanSignal){
	
MMM=dd$meanS
min=min(MMM)

for(i in 1:dim(MMM)[2]){
MMM[,i]=MMM[,i]+(abs(min)+offset)
} 

	MMM=log2(MMM)
	maintitle="MeanSignal"
	colorfill="orange"

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		dev.new()
		ddaux=dd
		ddaux$meanS=MMM	
		mvaMicroRna(ddaux,maintitle,verbose=FALSE)
		rm(ddaux)

		maintitle="MeanSignal- RLE "
		dev.new()
		RleMicroRna(MMM,maintitle,colorfill)

		dev.new()
		hierclusMicroRna(MMM,targets$GErep,methdis="euclidean",
        	methclu="complete",sel=FALSE,100)

}
# --- Gb="gProcessedSignal"  ----------

if(ProcessedSignal){
MMM=dd$procS
min=min(MMM)

for(i in 1:dim(MMM)[2]){
MMM[,i]=MMM[,i]+(abs(min)+offset)
} 
	MMM=log2(MMM) 
	maintitle="ProcessedSignal"
	colorfill="blue"
		
		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		dev.new()
		ddaux=dd
		ddaux$TPS=MMM	
		mvaMicroRna(ddaux,maintitle,verbose=FALSE)
		rm(ddaux)

		maintitle="ProcessedSignal - RLE "
		dev.new()
		RleMicroRna(MMM,maintitle,colorfill)

}
# --- Gf="gTotalProbeSignal"  ----------

if(TotalProbeSignal){
up=which(duplicated(dd$genes$ProbeName)==FALSE)
ddaux=dd[up,]
MMM=ddaux$TPS
min=min(MMM)

for(i in 1:dim(MMM)[2]){
MMM[,i]=MMM[,i]+(abs(min)+offset)
} 
	MMM=log2(MMM) 
	maintitle="TotalProbeSignal"
	colorfill="red"
		
		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		maintitle=" TotalProbeSignal - RLE "
		dev.new()
		RleMicroRna(MMM,maintitle,colorfill)

}
# --- Rf="gTotalGeneSignal"    ----------

if(TotalGeneSignal){
ddTGS=tgsMicroRna(dd,offset,half=FALSE,makePLOT=FALSE,verbose=FALSE)

	MMM=log2(ddTGS$TGS)
	maintitle="TotalGeneSignal"
	colorfill="green"

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)

		dev.new()
		plotDensityMicroRna(MMM,maintitle)

		maintitle=" TotalGeneSignal - RLE "
		dev.new()
		RleMicroRna(MMM,maintitle,colorfill)

}

# --- BGKmd="dd$other$gBGMedianSignal",     ----------

if(BGMedianSignal){

	MMM=log2(dd$other$gBGMedianSignal)
	maintitle="BGMedianSignal"
	colorfill="yellow"

		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)
}

# --- BGKus="gBGUsed")     ----------

if(BGUsed){

	MMM=log2(dd$other$gBGUsed)
	maintitle="BGused"
	colorfill="cyan"
		
		dev.new()
		boxplotMicroRna(MMM,maintitle,colorfill)
}

} # end function 

