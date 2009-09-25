\name{significantMicroRna}
\alias{significantMicroRna}
\title{ Summarize Differential Expression Results }
\description{
	The function summarizes the results from the differential expression analysis
	using the different objects that are obtained after 'limma' has been used, 
	such as the 'MArrayLM' object with the statistics and the 'TestResults' object
	highlighting the significant features.  
}
\usage{
significantMicroRna(eset, ddset, targets, fit2, 
CM, DE, DEmethod, MTestmethod, PVcut, Mcut,verbose=FALSE)
}
\arguments{
  \item{eset}{ ExpressionSet containing the Total Gene processed data }
  \item{ddset}{ An RGList object containing the Total Gene proceseed data }
  \item{targets}{data.frame with the target structure }
  \item{fit2}{MArrayLM object from \code{eBayes} 'limma' function}
  \item{CM}{Contrast matrix }
  \item{DE}{TestResults object}
  \item{DEmethod}{ method used in decideTests, only 'separate' or 'nestedF' are implemented }
  \item{MTestmethod}{method for multiple test }
  \item{PVcut}{p value threshold to declare significant features }
  \item{Mcut}{M value threshold to select within significant features}
  \item{verbose}{logical, if \code{TRUE} prints out output}
}
\details{
	A list containing the genes with their statistics is generated. The significant genes
	above the PVcut p values are also given in a html file that links the selected
	miRNAS to the miRBase  \url{http://microrna.sanger.ac.uk/}. A MA plots indicating the
	differentially expressed genes are also displayed. 
	
	When multiple contrasts are done, the method for the selection of the significant genes 
  	can be either 'separated' or 'nestedF'. See decideTests in
  	package limma \cite{limma} for a detailed description on these two methods. 
  	When 'separated' is used a list with all the genes that have been analized in 
  	limma is given. The list includes de following columns: 

		PROBE			- Probe name (one of the probes interrogating the gene)
		GENE    		- miRNA name
		PROBE chr\_coord 	- Agilent chromosomal location
		M       		- Fold change
		A       		- Mean of the intensity for that miRNA 
		t       		- moderated t-statistic
		pval    		- p value of the t-statistic
		adj.pval        	- p value adjusted by 'MTestmethod'
		fdr.pval		- p value adjusted by fdr

  Some times, the user can be set 'MTestmethod = none', in this case, it might be 
  interesting to still see the fdr value, despite of the fact that the user has 
  decided not apply any multiple testing correction. 

  If the 'nestedF' is used, then two lists are provided for each contrasts. A first
  containing the selected significant genes, and a second list containing the rest
  of the genes that have been analyzed. The columns given in this case is:

		PROBE			- Probe name (one of the probes interrogating the gene)
		GENE    		- miRNA name
		PROBE chr\_coord 	- Agilent chromosomal location
		M       		- Fold change
		A       		- Mean of the intensity for that miRNA 
		t       		- moderated t-statistic
		t pval    		- p value of the t-statistic
		F			- F statistic (null hypothesis: Ci = Cj, for all contrasts i, j)
		adj.F.pval        	- F p value adjusted by 'MTestmethod'
		fdr.F.pval		- F p value adjusted by fdr

  The html files, both for the 'separated' and 'nestedF' method, includes only the
  selected as significant genes. 
 

}
\references{ 
	Smyth, G. K. (2005). Limma: linear models for microarray data. In:
        'Bioinformatics and Computational Biology Solutions using R and
        Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.
        Huber (eds), Springer, New York, pages 397--420.
	
	miRBase: the home of microRNA data
	\url{http://microrna.sanger.ac.uk/}
 }
\author{ Pedro Lopez-Romero }
\seealso{
        An 'RGList' example containing proccesed data is in
        \code{ddPROC} and an overview of how the processed data is produced
        is given in \code{filterMicroRna}. The ExpressionSet object can be
        generated using \code{esetMicroRna}
	An overview of miRNA differential expression analysis is given in
        \code{basicLimma}
	An example of how to get the 'TestResults' object is in
        \code{getDecideTests}
}
\examples{
data(targets.micro)
data(ddPROC)
esetPROC=esetMicroRna(ddPROC,targets.micro,makePLOT=FALSE)

levels.treatment=levels(factor(targets.micro$Treatment))
treatment=factor(as.character(targets.micro$Treatment),
    levels=levels.treatment)
            
levels.subject=levels(factor(targets.micro$Subject))
subject=factor(as.character(targets.micro$Subject),
    levels=levels.subject)
	
design=model.matrix(~ -1 + treatment + subject  )

CM=cbind(MSC_AvsMSC_B=c(1,-1,0,0),
          MSC_AvsMSC_C=c(1,0,-1,0))
                
fit2=basicLimma(esetPROC,design,CM,verbose=TRUE)

DE=getDecideTests(fit2,
        DEmethod="separate",
        MTestmethod="BH",
        PVcut=0.10)

significantMicroRna(esetPROC,
	ddPROC,
	targets.micro,
  	fit2,
	CM,
	DE,
	DEmethod="separate",
	MTestmethod="BH",
	PVcut=0.10,
	Mcut=0,
	verbose=TRUE)
}
\keyword{documentation}
\keyword{utilities}

