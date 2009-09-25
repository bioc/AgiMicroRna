`paste.character` <-
function(paux){


         pi=paste(paux[1],paux[2],sep=",")
         if(length(paux) > 2 ){
         for(kk in 3:length(paux)){
         pi=paste(pi,paux[kk],sep=",")
	  } #kk
       } # if 
	return(pi)  
} # end function 

