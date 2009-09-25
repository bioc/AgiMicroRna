`readTargets` <-
function(infile,verbose=FALSE) {

	targets=read.table(file=infile,header=TRUE,
        	sep="\t",quote = "\"",dec=".",
        	strip.white=TRUE,blank.lines.skip=TRUE,
		fill=TRUE,skip=0,comment.char = "")

	if(names(targets)[1] != 'FileName'){
	stop('first column in target file must be called "FileName"',call. = FALSE)
	}
	if(names(targets)[2] != 'Treatment'){
	stop('second column in target file must be called "Treatment"',call. = FALSE)
	}
	if(names(targets)[3] != 'GErep'){
	stop('third column in target file must be called "GErep"',call. = FALSE)
	}

	rownames(targets)=targets$FileName  
	rownames(targets)=unlist(strsplit(rownames(targets),".txt"))
	if(verbose){
	cat("","\n")
	cat("Target File","\n")
	print(targets) 
	cat("","\n")
	}
	return(targets)

} # end function

