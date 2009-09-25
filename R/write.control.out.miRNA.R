`write.control.out.miRNA` <-
function(ddFILT,selSNR,targets){ 


	if (!is(ddFILT, "RGList")){
	  stop("'input' must be a RGList")
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
	values=round(ddFILT$G,3)
	flagsGID=matrix(ddFILT$other$gIsGeneDetected,nrow=dim(ddFILT)[1],ncol=dim(ddFILT)[2])

	GENE_ID=ddFILT$genes$GeneName
	probe.chr=ddFILT$other$chr_coord[,1]

		result=data.frame(as.character(PROBE_ID),as.character(GENE_ID),as.character(probe.chr),values)

		colnames(result)=c("PROBE","GENE","Probe Chr-Coord",paste(g2,g1,sep=" - "))

		outfile="NOCtrl_exprs.txt"
		write.table(result,file=outfile,row.names=F,
			col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")
		
		result=data.frame(PROBE_ID,as.character(GENE_ID),as.character(probe.chr),flagsGID)
		colnames(result)=c("PROBE & IsGeneDetected - (1 is Found)","GENE","Probe Chr-Coord",paste(g2,g1,sep=" - "))

		outfile="NOCtrl_FlagIsGeneDetected.txt"
		write.table(result,file=outfile,row.names=F,
			col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")

		rm(values)
} # end function 

