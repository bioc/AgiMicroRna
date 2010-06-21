`significantMicroRna` <-
function(eset,ddset,targets,fit2,
	CM,DE,DEmethod,MTestmethod,PVcut,Mcut,verbose=FALSE){

	require("geneplotter")
	require("codelink")  

	if (!is(eset, "ExpressionSet")){
	  stop("'eset' must be a ExpressionSet",call. = FALSE)
   	 	if (is.null(nrow(exprs(eset)))) {
        		stop("'eset' is empty",call. = FALSE)
	 	}
	}
	if (!is(ddset, "RGList")){
	  stop("'ddset' must be a RGList",call. = FALSE)
   	 	if (is.null(dim(ddset)[1])) {
        		stop("'ddset' is empty",call. = FALSE)
	 	}
	}
	if(missing(targets)){
		stop("'targets' is missing ",call. = FALSE)
	}else{
		if("GErep" %in% colnames(targets)){
			GErep=targets$GErep 
			nGE=sum(table(table(GErep)))
			g1=targets$GErep  # g1 must be numeric, from 1:n
			g2=rownames(targets)
		}else{
			stop("'targets' needs 'GErep' field",call. = FALSE)
		}
	}
	

	if (!is(fit2, "MArrayLM")){
	  stop("'fit2' must be a 'MArrayLM",call. = FALSE)
   	 	if (is.null(length(fit2$coefficients))) {
        		stop("fit2' is empty",call. = FALSE)
	 	}
	}
	if (!is(CM, "matrix")){
	  stop("'CM' must be a 'Contrast Matrix'",call. = FALSE)
   	 	if (is.null(dim(CM))) {
        		stop("'CM' is empty",call. = FALSE)
	 	}
	}
	if (!is(DE, "TestResults")){
	  stop("'DE' must be a 'TestResults'",call. = FALSE)
   	 	if (is.null(dim(DE))) {
        		stop("'DE' is empty",call. = FALSE)
	 	}
	}
	if(missing(DEmethod)){
		stop(" method for decideTests 'separate' or 'nestedF' is needed",call. = FALSE)
	}
	if(missing(MTestmethod)){
		stop(" method for multiple test 'none','BH' or 'BY' is needed",call. = FALSE)
	}
	if(missing(PVcut)){
		stop("'PVcut' is missing",call. = FALSE)
	} 	

	if(verbose){
	cat("------------------------------------------------------","\n")
	}
	
	nGEN=dim(eset)[1]
	nARR=dim(eset)[2]
	nCON=dim(CM)[2]


for(i in 1:nCON){ 

	if(verbose){
	cat("CONTRAST: ",i," - ",colnames(CM)[i],"\n")
	cat("\n")
	}

	notDE=which(DE[,i] == 0) # total number of DE genes
	nDDEE=which(DE[,i] != 0) # total number of DE genes 

   if(length(nDDEE) > 0){ 
	DEG.print.info(eset,DE,i,verbose)
   }
	
	
## SORTING objects by pvalue:

	# La seleccion por gen unico es por PVALOR. Se meten datos ordenados por PVALOR 
	method=match.arg(DEmethod,c("separate","nestedF"))
	switch(method,separate={

		adj.pval=p.adjust(fit2$p.value[,i],method=MTestmethod)
		fdr=p.adjust(fit2$p.value[,i],method="BH")
	
		# ordenamos file por pvalor 
		ord.pval=order(fit2$p.value[,i],decreasing=F)
		
		eset.ord=eset[ord.pval,]
		fit2.ord=fit2[ord.pval,]
		ddset.ord=ddset[ord.pval,]
		adj.pval=adj.pval[ord.pval]
		fdr=fdr[ord.pval]
	
	},nestedF={

		
		adj.Fpval=p.adjust(fit2$F.p.value,method=MTestmethod)
		fdr.Fpval=p.adjust(fit2$F.p.value,method="BH")
		
		# ordenamos file por F.pvalor:  separo nDDDEE de notDE 

		ord.pval.DE=nDDEE[order(fit2$F.p.value[nDDEE],decreasing=F)]
		ord.pval.notDE=notDE[order(fit2$F.p.value[notDE],decreasing=F)]

		eset.ord.DE=eset[ord.pval.DE,]
		fit2.ord.DE=fit2[ord.pval.DE,]
		ddset.ord.DE=ddset[ord.pval.DE,]
		adj.Fpval.DE=adj.Fpval[ord.pval.DE]
		fdr.Fpval.DE=fdr.Fpval[ord.pval.DE]

		eset.ord.notDE=eset[ord.pval.notDE,]
		fit2.ord.notDE=fit2[ord.pval.notDE,]
		ddset.ord.notDE=ddset[ord.pval.notDE,]
		adj.Fpval.notDE=adj.Fpval[ord.pval.notDE]
		fdr.Fpval.notDE=fdr.Fpval[ord.pval.notDE]

	})

	nDDEE.ord=seq(1:length(nDDEE))
	notDE.ord=seq((length(nDDEE)+1):dim(eset)[1])
		
## EXTRACTING RESULTS 

 	#	EXTRACTING M, A, T, AND PVALUES for the contrast [i]  - FOR THE WHOLE GENES 
	# 		EN SEPARATE SACO UNA UNICA LISTA CON TODOS LOS GENES, Y DE AQUI EL USER PUEDE
	#		ELEGIR SUS NDEE Y NOTDEE EN FUNCION DEL PVALOR. EN NESTEDF, LA SELECCION DE 
	#		GENES ES ALGO MAS COMPLEJA, Y SE SACAN LAS DOS LISTAS DE GENES, DEGS Y NOTDEGS 


method=match.arg(DEmethod,c("separate","nestedF"))

switch(method,separate={

## File Names 

	DEGLIST.ALL=paste("LIST.ALL",DEmethod,MTestmethod,
			colnames(CM)[i],"txt",sep=".")
## WRITE: ALL LIST.txt file: ALL genes 

	index.ALL=seq(1:dim(eset)[1])
	write.LIST.miRNA2(fit2.ord,index.ALL,DEGLIST.ALL,DEmethod,i,adj.pval,fdr,ddset.ord)

	filename=paste(DEGLIST.ALL,"html",sep=".")
 	genelist=fit2.ord$genes[adj.pval <= PVcut,]
 	title=paste("DEGs",DEmethod,MTestmethod,PVcut,colnames(CM)[i],sep=".")
 	head <- c("miRNA ID")
 	miRNA.htmlpage(genelist, filename, title, table.head=head,table.center = TRUE)

# MA PLOTS:
			MA.plot.miRNA(fit2,DE,CM,i)
	
},nestedF={

## File Names 
	
	DEGLIST=paste("DEGs",DEmethod,MTestmethod,PVcut,
			colnames(CM)[i],"txt",sep=".")
	NOTDE=paste("notDEG",Mcut,DEmethod,MTestmethod,PVcut,
			colnames(CM)[i],"txt",sep=".")

## WRITE: DEG.txt file: ALL nDDEE genes 
	# !!! PASAR LISTAS FIT ORDENADO, PARA QUE SAQUE LISTAS ORDENADAS POR PVALOR 
	
	write.LIST.miRNA2(fit2.ord.DE,nDDEE.ord,DEGLIST,DEmethod,i,adj.Fpval.DE,fdr.Fpval.DE,ddset.ord.DE)

	filename=paste(DEGLIST,"html",sep=".")
 	genelist=fit2.ord.DE$genes[nDDEE.ord,]
 	title=paste("DEGs",DEmethod,MTestmethod,PVcut,colnames(CM)[i],sep=".")
 	head <- c("miRNA ID")
 	miRNA.htmlpage(genelist, filename, title, table.head=head,table.center = TRUE) 

## WRITE: ALL notDEG.txt file: ALL notDEG genes 
# 
	write.LIST.miRNA(fit2.ord.notDE,notDE.ord,NOTDE,DEmethod,i,adj.Fpval.notDE,fdr.Fpval.notDE,ddset.ord.notDE)

## MA PLOTS:
		MA.plot.miRNA(fit2,DE,CM,i)

}) # END separate,nestedF


} # for(i in 1:nCON)
}  ## end extract.info.contrasts 

