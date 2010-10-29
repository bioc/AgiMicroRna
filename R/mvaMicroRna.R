`mvaMicroRna` <-
function(uRNAList,maintitle,verbose=FALSE) {


# MVA plots of each array against the genewise Median intensity

if(verbose){
cat("\n")
cat("------------------------------------------------------","\n")
cat("mvaMicroRna info:","\n")
}

nGEN=c(1:dim(uRNAList$meanS)[1])

indexrep=which(uRNAList$genes$ControlType == 0)

LL=length(indexrep)
if(verbose){
cat(" FEATURES :	",LL,"\n")
}
# ------------------ POS Ctrl --------------------
index1=which(uRNAList$genes$ControlType == 1)
index2=which(uRNAList$genes$GeneName[index1] != "DarkCorner")
index3=which(uRNAList$genes$GeneName[index1[index2]] != "miRNABrightCorner30")
index4=which(uRNAList$genes$GeneName[index1[index2[index3]]] != "SCorner3")

indexpos=nGEN[index1[index2[index3[index4]]]]

LL=length(indexpos)
if(verbose){
cat(" POSITIVE CTRL:		",LL,"\n")
}

## ------------------ NEG Ctrl --------------------
indexneg=which( uRNAList$genes$ControlType == -1)

LL=length(indexneg)
if(verbose){
cat(" NEGATIVE CTRL:		",LL,"\n")
}
## ------------------ STRUCTURAL  --------------------

indexDC=which(uRNAList$genes$GeneName == "DarkCorner" )
indexGE=which(uRNAList$genes$GeneName == "miRNABrightCorner30")
indexSC=which(uRNAList$genes$GeneName == "SCorner3")

LL=length(indexDC)+length(indexGE)+length(indexSC)
if(verbose){
cat(" STRUCTURAL:		",LL,"\n")
}

nARR=dim(uRNAList$meanS)[2]
y=apply(uRNAList$meanS,1,median)

for(i in 1:nARR) {

if (!missing(maintitle)){
what=paste(maintitle,colnames(uRNAList$meanS)[i]," - ","genewise Median")
}else{
what=paste(colnames(uRNAList$meanS)[i]," - ","genewise Median")
}

  x=uRNAList$meanS[,i]
		A=(x+y)/2
		M=(x-y)
	 mva=plot(A,M,type="p",cex=0.7,col="blue",xlab="A",ylab="M")

        points(A[indexrep],M[indexrep],cex=0.3,col="cyan3",pch=19)   	# Rep NonCtrl 
        points(A[indexDC],M[indexDC],cex=0.5,col="orange",pch=19)   	# DarkCorner  
        points(A[indexpos],M[indexpos],cex=0.5,col="red",pch=19)   	# pos Ctrl 
        points(A[indexneg],M[indexneg],cex=0.5,col="green",pch=19)   	# neg Ctrl 
        points(A[indexGE],M[indexGE],cex=1.0,col="yellow",pch=19)   	# miRNABrightCorner 
        points(A[indexSC],M[indexSC],cex=1.0,col="pink",pch=19)  		# SCorner3 

        smooth.fit = fitted(loess(M~A))
        points(A,smooth.fit,col='black',cex=0.5,pch=19)


        title(main=what)
        abline(0,0,col="black",lty=2, lwd=1)
        abline(2,0,col="navyblue")
        abline(-2,0,col="navyblue")

	colors=c("cyan","red","green","cyan3","orange","yellow","pink","black")
	samples=c("features","pos Ctrl","neg Ctrl","Rep NonCtrl","DarkCorner","miRNABrightCorner","SCorner3","M~A smooth fit")
	legend(x = "topright", legend = samples, cex = 0.8,fill = colors, inset = 0.05)
}

} # end function

