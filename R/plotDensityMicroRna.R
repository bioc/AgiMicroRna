`plotDensityMicroRna` <-
function(object,maintitle) {

	samples=colnames(object)
	nARR=dim(object)[2]

	colors <- rainbow(nARR,s=1,v=1,start = 0, 
		end = max(1,nARR - 1)/nARR, gamma = 1)
	y.max=c()
	x.max=c()

 	for (n in 1:nARR) {
 		y.max[n]=max(density(object[,n],na.rm=TRUE)$y)
 		x.max[n]=max(density(object[,n],na.rm=TRUE)$x)
 	}
 	y.pos=order(y.max,decreasing=TRUE,na.last=NA)
 	x.pos=order(x.max,decreasing=TRUE,na.last=NA)

 	for (n in y.pos) {
  		k=which(y.pos==n)
    		if (n ==y.pos[1])
         		plot(density(object[, n], na.rm = TRUE), 
			col = colors[n],main = "",
			asp=0.7*x.max[x.pos[1]]/y.max[y.pos[1]])
    		else lines(density(object[, n], na.rm = TRUE), col = colors[n])
 	}

	title(main=maintitle)
	legend(x = "topright", legend = samples, cex = 0.8,fill = colors, inset = 0.05)

} # end function

