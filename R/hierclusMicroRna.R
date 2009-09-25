`hierclusMicroRna` <-
function(object,GErep,methdis,methclu,sel,size){

samples=colnames(object) 


if(sel=="TRUE"){
  genes.var=apply(object,1,var)
  genes.var.select=order(genes.var,decreasing=TRUE)[1:size]
  object=object[genes.var.select,]
}

if(methdis=="euclidean"){
d=dist(t(object),method=methdis)  
main="HCLUST EUCLIDEAN"
}
if(methdis=="pearson"){
d=as.dist(1-cor(object,use="complete.obs",method=methdis)) # NA allowed 
main="HCLUST PEARSON"
}

if(sel=="TRUE"){
main=paste(main,"- high variance genes")
}else{
main=paste(main,"- all genes")
}

dim=dim(as.matrix(d))

 hc=hclust(d,method=methclu)
 plot(hc,labels=samples,main="")
 title(main=main)
}

