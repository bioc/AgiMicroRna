`mvaBASIC` <-
function(object,colorfill="blue",maintitle) {


nARR=dim(object)[2]
y=apply(object,1,median)

for(i in 1:nARR) {

if (!missing(maintitle)){
what=paste(maintitle,"-",colnames(object)[i]," vs. ","Median")
}else{
what=paste(colnames(object)[i]," vs. ","Median")
}

        x=object[,i]
        mva=plot((x+y)/2,(x-y),type="p",cex=.4,col=colorfill,xlab="A",ylab="M")
        title(main=what)
        abline(0,0,col="black",lty=2, lwd=1)
        abline(2,0)
        abline(-2,0)
}

} # end function 

