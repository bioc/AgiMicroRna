`write.LIST.miRNA` <-
function(fit2,index,outfile,DEmethod,i,adj.pval,fdr,ddset) {	

	if (!is(fit2, "MArrayLM")){
	  stop("'fit2' must be a ExpressionSet")	
	}
	
	if(!missing(ddset)){
		if(!is(ddset,"RGList")){
			stop("'ddset' must be a RGList")
		}
	}else{
		stop("'ddset' is needed")
	}

 method=match.arg(DEmethod,c("separate","nestedF"))
	switch(method,separate={

	GENE_ID=fit2$genes[index,]
	
	M=round(fit2$coefficients[index,i],3)
	A=round(fit2$Amean[index],3)
	t=round(fit2$t[index,i],3)
	pval=round(fit2$p.value[index,i],5)
	adj.pval=round(adj.pval[index],5)  
	fdr=round(fdr[index],5)            

		PROBE_ID=ddset$genes$ProbeName[index]

			result=data.frame(PROBE_ID,as.character(GENE_ID),
			M,A,t,pval,adj.pval,fdr)
			colnames(result)=c("PROBE","GENE",
			"M","A","t","pval","adj.pval","fdr.pval")
	 
	write.table(result,file=outfile,row.names=F,col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")

},nestedF={

	GENE_ID=fit2$genes[index,]
	
	M=round(fit2$coefficients[index,i],3)
	A=round(fit2$Amean[index],3)
	t=round(fit2$t[index,i],3)
	pval=round(fit2$p.value[index,i],5)

	F.stc=round(fit2$F[index],3)
	F.pval=round(fit2$F.p.value[index],5)
	adj.F.pval=round(adj.pval[index],5)  
	fdr.F.pval=round(fdr[index],5)

		PROBE_ID=ddset$genes$ProbeName[index]
	
			result=data.frame(PROBE_ID,as.character(GENE_ID),
			M,A,t,pval,F.stc,F.pval,adj.F.pval,fdr.F.pval)
		colnames(result)=c("PROBE","GENE",
		"M","A","t","t pval","F","F.pval","adj.F.pval","fdr.F.pval")
 
	write.table(result,file=outfile,row.names=F,col.names = TRUE,quote=F,dec=".",eol = "\n",sep = "\t")
})
} ## end function 

