`DEG.print.info` <-
function(eset,DE,i,verbose){
## DEGs
	
	if (!is(eset, "ExpressionSet")){
	  stop("'eset' must be a ExpressionSet")
   	 	if (is.null(nrow(exprs(eset)))) {
        		stop("'eset' is empty")
	 	}
	}
	if (!is(DE, "TestResults")){
	  stop("'DE' must be a 'TestResults'")
   	 	if (is.null(dim(DE))) {
        		stop("'DE' is empty")
	 	}
	}
	notDE=which(DE[,i] == 0) # total number of DE genes
	nDDEE=which(DE[,i] != 0) # total number of DE genes 
	posDE=which(DE[,i] > 0)
	negDE=which(DE[,i] < 0)

	cat("    DEGs: ",length(nDDEE),"\n")
	 
## POS DEGs: 
	
	cat("      UP: ",length(posDE),"\n")

## NEG DEGs: 

	cat("    DOWN: ",length(negDE),"\n")
	cat("\n")

} # end function

