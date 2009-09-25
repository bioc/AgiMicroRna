`getTDRows2` <-
function(ids){

# GIVES URLs in  "microrna.sanger.ac.uk" database 

# getQueryLink 
temp = paste("http://microrna.sanger.ac.uk/cgi-bin/sequences/mirna_entry.pl?acc=",ids,sep="")  

# getCells
out <- paste(" <A HREF=\"", temp, "\">", ids, "</A>",  sep = "")

# getTDRows 
rows <- paste("<TD>",out, "</TD>", sep = "")
return=(rows)
}

