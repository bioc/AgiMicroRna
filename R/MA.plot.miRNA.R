`MA.plot.miRNA` <-
function(fit2,DE,CM,i){
	posDE=which(DE[,i] > 0)
	negDE=which(DE[,i] < 0)

	M=fit2$coef[,i]
	A=fit2$Amean
      	dev.new()	

        plot(A,M,cex=1,col="black",pch=20)  # ALL OF THEM
	points(A,M,cex=0.6,col="cyan2",pch=19)  # ALL OF THEM

	points(A[posDE],M[posDE],cex=0.6,col="red",pch=19)    
	points(A[negDE],M[negDE],cex=0.6,col="blue",pch=19)    

        smooth.fit = fitted(loess(M~A))
        points(A,smooth.fit,col='black',cex=0.5,pch=19)

	colors=c("red","blue","black")
	samples=c("Up","Down","M~A smooth fit")
	legend(x = "topright", legend = samples, cex = 0.8,fill = colors, inset = 0.05)

	nn=paste(i," - ",colnames(CM)[i])
	title(main = nn)

} # end function 	

