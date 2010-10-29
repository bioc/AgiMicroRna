`write.control.out.miRNA` <-
function(ddFILT,selSNR,targets){ 


	if (!is(ddFILT, "uRNAList")){
	  stop("'input' must be a uRNAList")
   	 	if (is.null(dim(ddFILT)[1])) {
        		stop("'input' is empty")
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

	PROBE_ID=ddFILT$genes$ProbeName
	values=round(ddFILT$meanS,3)
	flagsGID=matrix(ddFILT$other$gIsGeneDetected,nrow=dim(ddFILT)[1],ncol=dim(ddFILT)[2])
	GENE_ID=ddFILT$genes$GeneName

			result=data.frame(as.character(PROBE_ID),as.character(GENE_ID),values)
			colnames(result)=c("PROBE","GENE",paste(g2,g1,sep=" - "))

		outfile="NOCtrl_exprs.txt"
		write.table(result,file=outfile,row.names=F,
			col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")
		
		result=data.frame(PROBE_ID,as.character(GENE_ID),flagsGID)
		colnames(result)=c("PROBE & IsGeneDetected - (1 is Found)","GENE",paste(g2,g1,sep=" - "))

		outfile="NOCtrl_FlagIsGeneDetected.txt"
		write.table(result,file=outfile,row.names=F,
			col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")

		rm(values)
} # end function 

