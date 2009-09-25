`miRNA.htmlpage` <-
function (genelist, filename, title, othernames, table.head, 
    table.center = TRUE) {

# modified from htmlpage, to use the getTDRows2 (modified to search in "microrna.sanger.ac.uk")


    outfile <- file(filename, "w")
    type <- "text/css"
    cat("<html>", "<head>", "<TITLE> miRNAs </TITLE>", 
        "</head>", "<body bgcolor=#FFFFFF >", "<H1 ALIGN=CENTER > miRNAs </H1>", 
        paste("<style type=", type, ">", sep = ""), "p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }", 
        "</style>", file = outfile, sep = "\n")
    if (!missing(title)) 
        cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n", 
            file = outfile, sep = "\n")
    if (table.center) 
        cat("<CENTER> \n", file = outfile)
    cat("<TABLE BORDER=4>", file = outfile, sep = "\n")
    if (!missing(table.head)) {
        headout <- paste("<TH>", table.head, "</TH>")
        cat("<TR>", headout, "</TR>", file = outfile, sep = "\n")
    }

  
    rows <- getTDRows2(genelist)
    nrows <- length(genelist)

    if (!missing(othernames)) {
        if (is.list(othernames)) {
            others <- ""
            for (nm in othernames) {
                if (is.matrix(nm)) {
                  for (i in 1:dim(nm)[2]) {
                    others <- paste(others, "<TD>", nm[, i], 
                      "</TD>", sep = "")
                  }
                }
                if (is.list(nm)) {
                  out <- vector()
                  for (j in seq(along = nm)) {
                    out[j] <- paste("<P>", nm[[j]], "</P>", sep = "", 
                      collapse = "")
                  }
                  out <- paste("<TD>", out, "</TD>", sep = "")
                  others <- paste(others, out, sep = "")
                }
                if ((is.vector(nm) || is.factor(nm)) && !is.list(nm)) 
                  others <- paste(others, "<TD>", nm, "</TD>", 
                    sep = "")
            }
        }
        else others <- paste("<TD>", othernames, "</TD>", sep = "")
        if (length(rows) != length(others)) 
            stop(paste("There are", length(rows), "rows in your genelist, but", 
                length(others), "rows in othernames.\n This will not give", 
                "good results!\n"))
        rows <- paste(rows, others)
    }
    for (i in 1:nrows) cat("<TR>", rows[i], "</TR>", file = outfile, 
        sep = "\n")
    cat("</TABLE>", file = outfile)
    if (table.center) 
        cat("</CENTER> \n", file = outfile)
    cat("</body>", "</html>", sep = "\n", file = outfile)
    close(outfile)
}

