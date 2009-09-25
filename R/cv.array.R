`cv.array` <-
function(ddDUP,rN,targets,what=c("ProbeName","GeneName"),verbose){


	what <- match.arg(what)

	Lreps=length(rN)
	nARR=dim(ddDUP)[2]

	cvmat=matrix(ncol=nARR,nrow=Lreps)
	colnames(cvmat)=rownames(targets) 
	rownames(cvmat)=rN

	if(what == "ProbeName"){	
			aux=ddDUP$genes$ProbeName
			maintitle=" CV of replicated Probes"
			colorfill="lightblue"
        		if(verbose){			
			cat("Replication at Probe level- MEDIAN  CV","\n")
			}
			
	}	
	if(what == "GeneName"){
			aux=ddDUP$genes$GeneName
			maintitle="  CV of replicated Genes "
			colorfill="lightgreen"
        		if(verbose){			
			cat("Replication at Gene level - MEDIAN  CV","\n")
			}
	}

	for(i in 1:Lreps){
  		index=which(aux==rN[i])
  		sd=round(apply(ddDUP$G[index,],2,sd),2)
  		m=round(apply(ddDUP$G[index,],2,mean),2)
  		cv=(sd/m)
  		cvmat[i,]=cv
 	}

	cvMedARR=round(apply(cvmat,2,median),3)

		boxplotMicroRna(cvmat,maintitle,colorfill,xlab="Samples",ylab="CV ")

        if(verbose){			
	print(cvMedARR)
	cat("------------------------------------------------------","\n")
	}

} # end function 

