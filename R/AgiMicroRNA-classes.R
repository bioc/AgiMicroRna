# AgiMicroRNA-classes.R defined from classes.R limma FILE (by Gordon Smyth)

setClass("uRNAList",
#  Class to hold Agilent MicroRna data
representation("list")
)

printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  from Gordon Smyth
{
	if(is.atomic(x)) {
		d <- dim(x)
		if(length(d)<2) which <- "OneD"
		if(length(d)==2) which <- "TwoD"
		if(length(d)>2) which <- "Array"
	} else {
		if(inherits(x,"data.frame")) {
			d <- dim(x)
			which <- "TwoD"
		} else {
			if(is.call(x))
				which <- "Call"
			else {
				if(is.recursive(x))
					which <- "Recursive"
				else
					which <- "Other"
			}
		}
	}
	switch(which,
	OneD={
		n <- length(x)
		if(n > 20) {
			print(x[1:5])
			cat(n-5,"more elements ...\n")
		} else
			print(x)
	},
	TwoD={
		n <- d[1]
		if(n > 10) {
			print(x[1:5,])
			cat(n-5,"more rows ...\n")
		} else
			print(x)
	},
	Array={
		n <- d[1]
		if(n > 10) {
			dn <- dimnames(x)
			dim(x) <- c(d[1],prod(d[-1]))
			x <- x[1:5,]
			dim(x) <- c(5,d[-1])
			if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
			dimnames(x) <- dn
			print(x)
			cat(n-5,"more rows ...\n")
		} else
			print(x)
	},
	Recursive={
		n <- length(x)
		if(n) {
			i <- names(x)
			if(is.null(i)) i <- seq_len(n)
			for (what in i) {
				y <- x[[what]]
				cat("$",what,"\n",sep="")
				Recall(y)
				cat("\n")
			}
		}
	},
	Call=,Other=print(x)
	)
}

setMethod("show","uRNAList",
#  Print and show method large data objects
#  Gordon Smyth
#  May 2003
function(object) {
        cat("An object of class \"",class(object),"\"\n",sep="")
        for (what in names(object)) {
                x <- object[[what]]
                cat("$",what,"\n",sep="")
                printHead(x)
                cat("\n")
        }
        for (what in setdiff(slotNames(object),".Data")) {
                x <- slot(object,what)
                if(length(x) > 0) {
                        cat("@",what,"\n",sep="")
                        printHead(x)
                        cat("\n")
                }
        }
})
## adapted from subsetting.R   (Gordon Smyth) # ~/R-packages/PACKAGES/EXAMPLES/limma/R
# [ method for subsetting the uRNALIST object 

setMethod("[", "uRNAList",
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

dim.uRNAList <- function(x) if(is.null(x$TGS)) c(0,0) else dim(as.matrix(x$TGS))
length.uRNAList <- function(x) prod(dim(x))
dimnames.uRNAList <- function(x) dimnames(x$TGS)

.setdimnames <- function(x, value)

{
	exprmatrices <- c("TGS","TPS","meanS","procS")
	for (a in exprmatrices) if(!is.null(x[[a]])) dimnames(x[[a]]) <- value
	for(a in names(x$other)) dimnames(x$other[[a]]) <- value
	if(!is.null(x$targets)) row.names(x$targets) <- value[[2]]
	x
}

"dimnames<-.uRNAList" <- .setdimnames
summary.uRNAList <- function(object,...) summary(unclass(object))



# onLoad stuff for S4 classes in NAMESPACE.
# .onLoad <- function(lib, pkg) {
#        require(methods, quietly = TRUE) || stop("Package methods unavailable!")

