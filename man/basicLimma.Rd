\name{basicLimma}
\alias{basicLimma}
\title{ Linear models Using limma  }
\description{
  Differential expression analysis using the linear model features implemented
  in the limma package. A linear model is fitted to each miRNA gene so that the 
  fold change between different experimental conditions and their standard 
  errors can be estimated. Empirical Bayes methods are applied to obtain moderated
  statistics 
}
\usage{
basicLimma(eset, design, CM,verbose = FALSE)
}
\arguments{
  \item{eset}{ExpressionSet containing the processed log-expression values }
  \item{design}{design matrix }
  \item{CM}{contrast matrix }
  \item{verbose}{logical, if \code{TRUE} prints out output}
}
\details{

  In our data example (see the target file in Table 1 in vignette), 
  we have used a paired design (by subject) to assess the differential expression 
  between two treatments B and C vs a control treatment A. That is, we want to obtain the 
  microRNAS that are differentially expressed between conditions A vs B and A vs C. 
  The linear model that we are going to fit to every miRNA is defined by equation:  
  y = Treatment + Subject + error term. This model is going to 
  estimate the treatment effect and  then, the comparison between the different treatments 
  are done in terms of contrasts between the estimates of the treatment effects. 
  To fit the model, we  need first to define a design matrix. The design matrix is an 
  incidence matrix that relates each array/sample/file to its given experimental
  conditions, in our case, relates each file to one of the three treatments and 
  with its particular subject. 
  If treatment is a factor variable, we can define de desing matrix using 
  model.matrix(~ -1 + treatment + subject). Then the linear model can be fitted using 
  fit=lmFit(eset,design). This will get the treatment estimates for each 
  microRNA in the eset object: 
  
               treatmentA treatmentB treatmentC subject2
hsa-miR-152        7.5721      7.656      7.566  -0.1157
hsa-miR-15a*       0.9265      1.066      1.211  -0.2242
hsa-miR-337-5p     6.2448      7.298      7.084  -0.4489

	
 We can define the contrasts  of interest using a contrast matrix as in 
		CM=cbind(MSC\_AvsMSC\_B=c(1,-1,0),
                	MSC\_AvsMSC\_C=c(1,0,-1))

 And then, we can estimate those contrats using fit2=contrasts.fit(fit,CM). Finally, we
 can obtain moderated statistics using fit2=eBayes(fit2). 
	
 The function 'basicLimma' implemented in AgiMicroRna produces the last fit2 
 object, that has in fit2\$coeff the M values, in fit\$t the moderated-t 
 statistic of the contrasts, and in fit2\$p.value the corresponding
 p value of each particular contrasts. Be aware that these p values must be corrected by
 multiple testing.  

	               MSC\_AvsMSC\_B MSC\_AvsMSC\_C
	hsa-miR-152      0.67567761  0.977326746
	hsa-miR-15a*     0.68019442  0.413657270
	hsa-miR-337-5p   0.03737814  0.075248741	

 See limmaUsersGuide() for a complete description of the limma package. 

}
\value{
	An MArrayLM object of the package limma
}
\references{ 
	Smyth, G. K. (2005). Limma: linear models for microarray data. In:
  	'Bioinformatics and Computational Biology Solutions using R and
  	Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.
  	Huber (eds), Springer, New York, pages 397--420.
	
	Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing
	diferential expression in microarray experiments. Statistical Applications in 
	Genetics and Molecular Biology, Vol. 3, No. 1, Article 3.
	http://www.bepress.com/sagmb/vol3/iss1/art3	

}
\author{ Pedro Lopez-Romero }
\seealso{
        An 'RGList' example containing proccesed data is in
        \code{ddPROC} and an overview of how the processed data is produced
        is given in \code{filterMicroRna}. The ExpressionSet object can be
        generated using \code{esetMicroRna}
}
\examples{
\dontrun{
data(targets.micro)
data(ddPROC)
esetPROC=esetMicroRna(ddPROC,targets.micro,makePLOT=FALSE,verbose=FALSE)

levels.treatment=levels(factor(targets.micro$Treatment))
treatment=factor(as.character(targets.micro$Treatment),
  		levels=levels.treatment)
            
levels.subject=levels(factor(targets.micro$Subject))
subject=factor(as.character(targets.micro$Subject),
  		levels=levels.subject)
	
design=model.matrix(~ -1 + treatment + subject   )

CM=cbind(MSC_AvsMSC_B=c(1,-1,0,0),
          MSC_AvsMSC_C=c(1,0,-1,0))
                
fit2=basicLimma(esetPROC,design,CM,verbose=TRUE)

names(fit2)
head(fit2$coeff)
head(fit2$p.value)
plot(fit2\$Amean,fit2$coeff[,1],xlab="A",ylab="M")								
abline(h=0)
abline(h=c(-1,1),col="red")
plot(fit2$coeff[,1],fit2$p.value[,1], xlab="M",ylab="p value")
}

}
\keyword{documentation}
\keyword{utilities}
