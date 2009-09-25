`boxplotMicroRna` <-
function(object,maintitle="",colorfill="blue",xlab="samples",ylab="expression") {

opar <- par(las = 3) 
    on.exit(par(opar))

	nARR=dim(object)[2]
	med=apply(object,2,mean)

        MIN=min(object,na.rm=TRUE)
        MAX=max(object,na.rm=TRUE)

	plot(med,xlim=c(0,nARR+1),ylim=c(MIN,MAX),
	axes=FALSE,xlab=xlab,ylab=ylab)

        boxplot(data.frame(object),col=colorfill,outline=TRUE,add=TRUE)
        points(med,type = "p",col ="blue")
        lines(med,type = "l",col ="blue",lty = "dotted")
	title(main=maintitle)

} # end function

