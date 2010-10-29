`rmaMicroRna` <-

function(dd,normalize, background){
    
require(affy)
require(preprocessCore)

if (!is(dd, "uRNAList")){
     stop("'input' must be a uRNAList")
          if (is.null(dim(dd)[1])) {
                   stop("'input' is empty")
          }
}

# expression matrix  "gMeanSignal = dd$meanS"
yy= dd$meanS

# background correction
if(background == TRUE ){

yy=rma.background.correct(yy,copy=TRUE)

  min=min(yy)
    for(i in 1:dim(yy)[2]){
      yy[,i]=yy[,i]+(abs(min)+ 2)
    }  
}

# quantiles
if(normalize == TRUE){
	dd$meanS=normalizeBetweenArrays(yy,method='quantile')
}else{ 
	dd$meanS = yy
}

# probe summarization / median
ProbeName.rep = unique(dd$genes$ProbeName[duplicated(dd$genes$ProbeName)])
dd.aux=dd[1:length(ProbeName.rep),]

for(ii in 1: length(ProbeName.rep)){
 index=which(dd$genes$ProbeName %in% ProbeName.rep[ii])
 dd.aux$meanS[ii,]=apply(dd$meanS[index,],2,median)

 dd.aux$genes[ii,] = dd$genes[index[1],]
 dd.aux$other$gIsGeneDetected[ii,] = dd$other$gIsGeneDetected[index[1],]
 dd.aux$other$gIsSaturated[ii,] = dd$other$gIsSaturated[index[1],]
 dd.aux$other$gIsFeatNonUnifOL[ii,] = dd$other$gIsFeatNonUnifOL[index[1],]
 dd.aux$other$gIsFeatPopnOL[ii,] = dd$other$gIsFeatPopnOL[index[1],]
 dd.aux$other$chr_coord[ii] = dd$other$chr_coord[index[1]]
}
 
if(min(dd.aux$meanS) < 0){
dd.aux$meanS = dd.aux$meanS + abs(min(dd.aux$meanS))+ 0.5
}

# rma 
pNList = dd.aux$genes$GeneName   # character with (replicated) Genes
ngenes <- length(unique(pNList))
pNList <- split(0:(length(pNList) - 1), pNList)  	
				# list with names of probes, and positions in yy/pNList

exprs <- .Call("rma_c_complete_copy", dd.aux$meanS, pNList, ngenes, 
              normalize=FALSE, background=FALSE,bgversion=2,
              verbose=TRUE, PACKAGE = "affy")

ddTGS.rma=dd.aux[1:length(rownames(exprs)),] 

for(ii in 1: length(rownames(exprs))){

 index=which(dd.aux$genes$GeneName %in% rownames(exprs)[ii])
  ddTGS.rma$meanS[ii,]=exprs[ii,]

 ddTGS.rma$genes[ii,] = dd.aux$genes[index[1],]
 ddTGS.rma$other$gIsGeneDetected[ii,] = dd.aux$other$gIsGeneDetected[index[1],]
 ddTGS.rma$other$gIsSaturated[ii,] = dd.aux$other$gIsSaturated[index[1],]
 ddTGS.rma$other$gIsFeatNonUnifOL[ii,] = dd.aux$other$gIsFeatNonUnifOL[index[1],]
 ddTGS.rma$other$gIsFeatPopnOL[ii,] = dd.aux$other$gIsFeatPopnOL[index[1],]
 ddTGS.rma$other$chr_coord[ii] = dd.aux$other$chr_coord[index[1]]
}

ddTGS.rma$TGS = ddTGS.rma$meanS
ddTGS.rma$TGS = ddTGS.rma$meanS
ddTGS.rma$TGS = ddTGS.rma$meanS 
return(ddTGS.rma)
}
