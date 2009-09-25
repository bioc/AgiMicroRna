`getDecideTests` <-
function(fit2,DEmethod,MTestmethod,PVcut,verbose=FALSE){


	if (!is(fit2, "MArrayLM")){
	  stop("'design' must be a 'MArrayLM")
   	 	if (is.null(length(fit2$coefficients))) {
        		stop("fit2' is empty")
	 	}
	}

	if(missing(DEmethod)){
		stop(" method for decideTests 'separate' or 'nestedF' is needed")
	}
	if(missing(MTestmethod)){
		stop(" method for multiple test 'none','BH' or 'BY' is needed")
	}
	if(missing(PVcut)){
		stop("'PVcut' is missing")
	} 	


	DE=decideTests(fit2,method=DEmethod
		,adjust.method=MTestmethod,p.value=PVcut)

	sumDE= summary(DE)
	rownames(sumDE)[1]="DOWN"
	rownames(sumDE)[3]="UP"
	sumDE=sumDE[c(3,1),]

	if(verbose){
	cat("\n")
	cat("------------------------------------------------------","\n")
	cat(" Method for Selecting DEGs:",DEmethod,"\n")
	cat(" Multiple Testing  method: ",MTestmethod,"- pval",PVcut,"\n")
	# cat(" DEG Significants:",length(ord),"\n") 
      	cat("\n")
      	print(sumDE)
      	cat("------------------------------------------------------","\n")
	}

	return(DE)
 } # end function 

