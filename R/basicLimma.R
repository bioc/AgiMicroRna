`basicLimma` <-
function(eset,design,CM,verbose = FALSE) { 

	if (!is(eset, "ExpressionSet")){
	  stop("'eset' must be a ExpressionSet", call. = FALSE)
   	 	if (is.null(nrow(exprs(eset)))) {
        		stop("'eset' is empty")
	 	}
	}

	if (!is(design, "matrix")){
	  stop("'design' must be a 'Design Matrix'",call. = FALSE)
   	 	if (is.null(dim(design))) {
        		stop("design' is empty",call. = FALSE)
	 	}
	}

	if (!is(CM, "matrix")){
	  stop("'CM' must be a 'Contrast Matrix'",call. = FALSE)
   	 	if (is.null(dim(CM))) {
        		stop("'CM' is empty",call. = FALSE)
	 	}
	}

if(verbose){
	cat("DATA","\n")
	print(dim(eset))
	cat("\n")
 }
	fit=lmFit(eset,design) 
	fit2=contrasts.fit(fit,CM)
	fit2=eBayes(fit2)

	return(fit2)

} # end function 

