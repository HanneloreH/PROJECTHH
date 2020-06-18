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
library(optrees)
library(ape)

### GET DATA
setwd("~/Shared Folder Hannelore/CENTOS PROJECT HH SHARE/Making Tree in R")
cgtable <- read.csv(file = "results_allelesKP_2020test.tsv",
                        sep = "\t", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
metadata <- read.csv(file = "MetadataKP_2020test.csv",
                     header = TRUE, 
                     stringsAsFactors = FALSE)
### GETTING TO KNOW THE DATA
dim(cgtable)
dim(metadata)
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
#???https://sites.google.com/site/therepiproject/r-pac/about



#installing outbreaktools
#problems!!





################################################# ALTERNATIVE 5: Incidence/outbreaker/...

g <- transGraph(cgtable, thres=0)



################################################# ALTERNATIVE 6: GrapeTRee



################################################# ALTERNATIVE 7: getMinimumSpanningTree
#http://www.smartana.co.uk/mst-tutorial/MST.html


dim(cgtable)
dim(metadata)

library(optrees)
library(igraph)
library(ape)
library(networkD3)

#make connected weighted undirected graph to get arcs:
#arcs = matrix with the list of arcs of the graph. Each row represents one arc. The first two columns contain the two endpoints 
#of each arc and the third column contains their weights.
#change rownames
row.names(cgtable)<-cgtable$FILE
#remove first column
cgtable<-cgtable[,-1]
cgtable[,1:3]
#make matrix
cgmatrix <- as.matrix(cgtable)
######make adjacency matrix
cgmatrix[,1:3]
#distance matrix: pairwise (alleles)
cgdist <-dist.gene(cgmatrix, method= "pairwise", pairwise.deletion=FALSE)
cgdist
#remark first sample is not in rows but in columns (probably because otherwise empty row?)
#terug in matrix omzetten
cgdist <-as.matrix(cgdist)
#make graph
cggraph = graph_from_adjacency_matrix(cgdist, mode="undirected", weighted=TRUE, diag=TRUE)
plot(cggraph)
cgmst = mst(cggraph) # no longer works??
plot(cgmst)
cgmst2 = minimum.spanning.tree(cggraph) #does work, not the same?
plot(cgmst2)
#add weights
plot(cgmst2, weights=edgematrix )


#PLOT LAYOUT
#https://kateto.net/network-visualization
install.packages("igraph") 
install.packages("network") 
install.packages("sna")
install.packages("ggraph")
install.packages("visNetwork")
install.packages("threejs")
install.packages("networkD3")
install.packages("ndtv")

library("igraph") 
library("network") 
library("sna")
library("ggraph")
library("visNetwork")
library("threejs")
library("networkD3")
library("ndtv")

#nodes and links/edges
cgnetwork= igraph_to_networkD3(cggraph)
nodes = cgnetwork$nodes #hier is staal SA32 = 1
links = cgnetwork$links #hier is staal SA32 = 0
linksValue = links$value
#igraph object
net <- graph_from_adjacency_matrix(cgdist, mode="undirected", weighted=TRUE, diag=TRUE)
net
#UNW = undirected, named, weighted, 16 nodes, 120 edges
#see edges
E(net)
#see nodes
V(net)
net[5,7]
#get edgelist
edgelist = as_edgelist(net,names=T)
edgematrix = as_adjacency_matrix(net, attr="weight")
#simplify plot
plot(net)
net2= simplify(net, remove.multiple=F, remove.loops=T)
plot(net, layout=layout_with_mds)
plot(net2, layout=layout_with_mds)



##############################################################################
#tree with ALL links

#transform into readable format for networkD3
cgnetwork= igraph_to_networkD3(cggraph)
#interessant een deel links: met alle arcs
#en een deel nodes: met alle nodes
#plot
simpleNetwork(cgnetwork$links, zoom=TRUE)


###################################################################################
# Classical MDS
# N rows (objects) x p columns (variables)
# distance table
cgdist
#make a fit(?)
fit <- cmdscale(cgdist,eig=TRUE, k=2) # k is the number of dim
# plot solution 
plot(fit$points[,1], fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", 
     main="KP-data MDS",    type="n")
text(x, y, labels = row.names(cgtable), cex=.7)
#focus on outbreak
plot(fit$points[,1], fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", 
     main="KP-data MDS DETAIL",    type="n", ylim=c(-10,-4), xlim=c(-530,-515) )
text(x, y, labels = row.names(cgtable), cex=.7)
###################################################################################





