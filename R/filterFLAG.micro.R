`filterFLAG.micro` <-
function(flag,GErep,nGE,minFLAG,limSNR) {
    
        index=0
        ii=0
        while(ii < nGE) {
        ii=ii+1
        minSize=length(which(GErep == unique(GErep)[ii]))*(limSNR/100)  
        pos=which(GErep == ii)  
                                
        aux=flag[pos]
        nan=!is.na(aux)         
        aux=aux[nan]           
        lp=length(aux)         

                sumSNR=0
                jj=0
                while(jj < lp ){
                jj=jj+1

                        if(aux[jj] == minFLAG){
                        sumSNR=sumSNR+1
                                if(sumSNR >= minSize){
                                index=1
                                jj=lp + 1  # exit the while(jj)
                                ii=nGE +1  # exit the while(ii)
                                }
                        }
                }
         }
        
	return(index)
} # end of function 

