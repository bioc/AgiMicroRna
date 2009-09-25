`writeEset` <-
function(eset,ddPROC,targets,verbose=FALSE) {


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
		GErep=targets$GErep 
		nGE=sum(table(table(GErep)))
		g1=targets$GErep  # g1 must be numeric, from 1:n
		g2=rownames(targets)
	}else{
		stop("'targets' needs 'GErep' field")
	}
	
	GENE_ID=featureNames(eset)
	values=round(exprs(eset),3)
	PROBE_ID=ddPROC$genes$ProbeName
	chr_coord=ddPROC$other$chr_coord[,1]
	

		result=data.frame(PROBE_ID,as.character(GENE_ID),
			as.character(chr_coord),values)

		colnames(result)=c("PROBE","GENE","Probe chr_coord",paste(g2,g1,sep=" - "))

	outfile="ProcessedData.txt"
	write.table(result,file=outfile,row.names=F,col.names = TRUE,
		quote=F,dec=".",eol = "\n",sep = "\t")
	
	if(verbose){
	cat("PROCESSED DATA  :",length(PROBE_ID),"\n")
	cat("\n")
	}
} # END.FUNCTION 

