## adapted from subsetting.R   (Gordon Smyth) # ~/R-packages/PACKAGES/EXAMPLES/limma/R
# [ method for subsetting the uRNALIST object 

assign("[.uRNAList",
# setMethod("[", "uRNAList",
function(x, i, j, ..., drop=FALSE) {

        if(!missing(i)) {
                x$TGS = x$TGS[i,]
                x$TPS = x$TPS[i,]
                x$meanS = x$meanS[i,]
                x$procS = x$procS[i,]

                # genes (1d) $genes is data.frame 
                x$genes = x$genes[i,]

                # other (2d) 
                x$other$gIsGeneDetected = x$other$gIsGeneDetected[i,]
                x$other$gIsSaturated = x$other$gIsSaturated[i,]
                x$other$gIsFeatNonUnifOL = x$other$gIsFeatNonUnifOL[i,]
                x$other$gIsFeatPopnOL = x$other$gIsFeatPopnOL[i,]
                x$other$gBGMedianSignal = x$other$gBGMedianSignal[i,]
                x$other$gBGUsed = x$other$gBGUsed[i,]
        }
        if(!missing(j)) {
                x$TGS <- x$TGS[,j]
                x$TPS <- x$TPS[,j]
                x$meanS <- x$meanS[,j]
                x$procS <- x$procS[,j]

                # other (2d) 
                x$other$gIsGeneDetected = x$other$gIsGeneDetected[,j]
                x$other$gIsSaturated = x$other$gIsSaturated[,j]
                x$other$gIsFeatNonUnifOL = x$other$gIsFeatNonUnifOL[,j]
                x$other$gIsFeatPopnOL = x$other$gIsFeatPopnOL[,j]
                x$other$gBGMedianSignal = x$other$gBGMedianSignal[,j]
                x$other$gBGUsed = x$gBGUsed[,j]

                x$targets = x$targets[j,1]
        }
        return(x)
})
