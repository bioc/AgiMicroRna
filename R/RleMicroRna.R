`RleMicroRna` <-
function(object,maintitle="",colorfill="blue") {


nARR=dim(object)[2]
nGEN=dim(object)[1]

y=apply(object,1,median)
mva=matrix(nrow=nGEN,ncol=nARR)

for(i in 1:nARR) {
        x=object[,i]
        mva[,i]=(x-y)
}

	med=apply(mva,2,median)
        MIN=min(mva,na.rm=TRUE)
        MAX=max(mva,na.rm=TRUE)

        par(las=3)
	plot(med,xlim=c(0,nARR+1),ylim=c(MIN,MAX),
	axes=FALSE,xlab="Samples",ylab="M")

	colnames(mva)=colnames(object)
	boxplot(data.frame(mva),outline=TRUE,add=TRUE,col=colorfill)

	points(med,type = "p",col = "blue")
	lines(med,type = "l",col = "blue",lty = "dotted")
	title(main=maintitle)
	
	abline(0,0,col="red")
	par(las=0)

} # end function 

