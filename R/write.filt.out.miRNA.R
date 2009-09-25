`write.filt.out.miRNA` <-
function(ddFILT,selSNR,outfile,FLAG,targets){ 


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

	PROBE_ID=ddFILT$genes$ProbeName[-selSNR]  
	FLAG=FLAG[-selSNR,]	
	GENE_ID=ddFILT$genes$GeneName[-selSNR]
	probe.chr=ddFILT$other$chr_coord[-selSNR,1]

	if(length(PROBE_ID) != 0){
 	if(length(PROBE_ID) ==1){

 		result=data.frame(PROBE_ID,as.character(GENE_ID),
			as.character(probe.chr),t(FLAG))
	}
 	if(length(PROBE_ID) > 1){
 		result=data.frame(PROBE_ID,as.character(GENE_ID),
			as.character(probe.chr),FLAG)
	}
		colnames(result)=c("PROBE","GENE","Probe Chr-Coord",paste(g2,g1,sep=" - "))

		write.table(result,file=outfile,row.names=F,
			col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")
	}
		rm(FLAG)
} # end function 

