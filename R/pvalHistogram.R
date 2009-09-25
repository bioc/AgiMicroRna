`pvalHistogram` <-
function(fit2,DE,PVcut,DEmethod,MTestmethod,CM,verbose=FALSE){

	nCON=dim(DE)[2]

	if(nCON > 1 && DEmethod == "nestedF"){
	## selecciona F.p.value corrigiendo por test-multiple 

		ord=which(p.adjust(fit2$F.p.value,method=MTestmethod) <= PVcut)
		if(verbose){
		cat("num F Contrasts:",nCON,"\n")
		cat("DEG Significants by F.p.value <= ",PVcut,"-",length(ord),MTestmethod,"\n") 
		cat("------------------------------------------------------","\n")
		}
		aux=paste.character(colnames(CM))
		maintitle=paste("F-test ",aux,DEmethod,sep=" - ")
		hist(fit2$F.p.value,freq=TRUE,col="yellow",
			main=maintitle,xlab="F test p.value",ylab=" freq ")

	}
	if(nCON == 1 || DEmethod == "separate"){

		for(i in 1:nCON){ 
			maintitle=paste(colnames(CM)[i],DEmethod,sep=" - ")
			dev.new()
			hist(fit2$p.value[,i],freq=TRUE,col="green",
				main=maintitle,xlab="t test p.value",ylab=" freq ")
		}
	}
} # end function 

