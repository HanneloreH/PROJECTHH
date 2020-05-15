#Making a minimal spanning tree in R



### LOAD REQUIRED PACKAGES/LIBRARIES
library(igraph) #for MST
library(SIT) #for MST
library(limma)
library(edgeR) #eigenlijk voor RNA-seq expression profiles (maar is gelijkaardig :) )
library(adegenet) #om genind object te maken
library(OutbreakTools) #no longer on CRAN because not maintained, but downloaded most recent from archive
library(incidence)
library(outbreaker)
library(poppr)
library(network)


### GET DATA
setwd("~/Shared Folder Hannelore/CENTOS PROJECT HH SHARE")
cgtable <- read.csv(file = "cgMLST.tsv",
                        sep = "\t", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
metadata <- read.csv(file = "metadata.csv",
                     header = TRUE, 
                     stringsAsFactors = FALSE)
### GETTING TO KNOW THE DATA
dim(cgtable)
class(cgtable)
class(metadata)




################################################# ALTERNATIVE 1: igraph

#installing
#edgeR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

## MAKE DISTANCE MATRIX -> FOUT: neemt in rekening het nummer van de type van de loci, maar is niet nodig, moet gewoon verschil aanduiden...
cgDist2<-as.matrix(dist(c
cgGraph <- graph_from_data_frame(cgDist2, directed=TRUE, vertices=NULL)gtable))
dim(cgDist2)
#creating an igraph graph from a data frame
class (cgGraph)
#do mst
mst(cgGraph, weights=NULL, algorithm=NULL)
#only text??


################################################ ALTERNATIVE 2 : SIT
#INSTALL
library(curl)
library(quantmod)
library(devtools)
library(zoo)
devtools::install_github('systematicinvestor/SIT.date')
curl_download('https://github.com/systematicinvestor/SIT/raw/master/SIT.tar.gz', 'sit',mode = 'wb',quiet=T)
install.packages('sit', repos = NULL, type='source')

library(SIT)
load.packages("quantmod")

ret = diff(log(cgtable))

tickers = nasdaq.100.components() #ERROR




################################################# ALTERNATIVE 3: 
#NO: because specific on genotypes
# https://grunwaldlab.github.io/poppr/reference/poppr.msn.html 


#FUNCTION diss.dist: Calculate a distance matrix based on relative dissimilarity
diss.dist <- function(x, percent=FALSE, mat=FALSE){
  stopifnot(is(x, "gen"))
  ploid     <- x@ploidy
  if (is(x, "bootgen")){
    ind.names <- x@names
  } else {
    ind.names <- indNames(x)
  }
  inds      <- nrow(x@tab)
  np        <- choose(inds, 2)
  dist.mat  <- matrix(data = 0, nrow = inds, ncol = inds)
  numLoci   <- nLoc(x)
  type      <- x@type
  if (type == "PA"){
    dist_by_locus <- matrix(.Call("pairdiffs", x@tab))
    ploid <- 1
  } else if (is(x, "bootgen")){
    dist_by_locus <- vapply(seq(numLoci), function(i){
      .Call("pairdiffs", tab(x[, i]))/2
    }, numeric(np))
  } else {  
    x <- seploc(x)
    dist_by_locus <- vapply(x, function(x) .Call("pairdiffs", x@tab)/2,
                            numeric(np))
  }
  if (is.matrix(dist_by_locus)){
    dist.mat[lower.tri(dist.mat)] <- rowSums(ceiling(dist_by_locus))    
  } else {
    dist.mat[lower.tri(dist.mat)] <- ceiling(dist_by_locus)
  }
  colnames(dist.mat) <- ind.names
  rownames(dist.mat) <- ind.names
  if (percent){
    dist.mat <- sweep(dist.mat, 1, ploid * numLoci, "/")
  }
  dist.mat <- as.dist(dist.mat)
  if (mat == TRUE){
    dist.mat <- as.matrix(dist.mat)
  }
  return(dist.mat)
}


#make genind object
cgtableG <- is.genind(cgtable)
cgDist3 <- diss.dist(cgtableG)



################################################# ALTERNATIVE 4:  https://rdrr.io/cran/OutbreakTools/man/plotggMST.html
#♦https://sites.google.com/site/therepiproject/r-pac/about



#installing outbreaktools
#problems!!





################################################# ALTERNATIVE 5: Incidence/outbreaker/...

g <- transGraph(cgtable, thres=0)



################################################# ALTERNATIVE 5: GrapeTRee




