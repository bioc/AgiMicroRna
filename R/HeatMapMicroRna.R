`HeatMapMicroRna` <-
function(object,size,maintitle) {

require(marray) # maPalette
require(gplots) # heatmap.2 
require(gtools)
require(gdata)

# ...........................................
samples=colnames(object)
names=rownames(object)
# ...........................................
genes.var=apply(object,1,var,na.rm=TRUE)
genes.var.select=order(genes.var,decreasing=TRUE)[1:size]
DD.s= object[genes.var.select,]

samples=colnames(DD.s)
names=rownames(DD.s)

c <- rainbow(ncol(DD.s), start=0, end=.3)
rc <- rainbow(nrow(DD.s), start=0, end=.3)

# col=gentlecol(256)
#col <- colorRampPalette(c("green", "red"))(123)
rbg=maPalette(low="green",high="red",mid="black",k=50)

heatmap.2(DD.s,labCol=samples,labRow=names,scale="none",
	col=rbg,margin=c(10,10),tracecol="cyan")

		#scale=row : normalize by rows
		#scale=col : normalize by cols 

subtitle=paste(as.character(size)," high variance genes")
        title(main =maintitle,sub=subtitle)
# ...........................................

} # end function

