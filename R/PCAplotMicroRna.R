`PCAplotMicroRna` <-
function(eset,targets) {

require(affycoretools)

	if (!is(eset, "ExpressionSet")){
	  stop("'object' must be a ExpressionSet")
   	 	if (is.null(nrow(exprs(eset)))) {
        		stop("'object' is empty")
	 	}
	}
	if(missing(targets)){
		stop("'targets' is missing ")
	}

	if("GErep" %in% colnames(targets)){
		g1=targets$GErep  # g1 must be numeric, from 1:n
		g2=rownames(targets)
	}else{
		stop("'targets' needs 'GErep' field")
	}
	DF=data.frame(exprs(eset))
	m=as.matrix(na.omit(DF))
	exprs(eset)=m

	plotPCA(eset,g1,g2,addtext=g2,screeplot=FALSE,
		x.coord=NULL,y.coord=NULL,squarepca=FALSE)

	dev.new()
	plotPCA(eset,g1,g2,screeplot=TRUE,x.coord=NULL,y.coord=NULL)

rm(eset)
rm(DF)
rm(m)

} # end function
