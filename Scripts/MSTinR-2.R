
### LOAD REQUIRED PACKAGES/LIBRARIES
library(optrees)
library(igraph)
library(ape)
library(networkD3)
library(RColorBrewer)
library(ggplot2)
library(emstreeR)


### GET DATA
setwd("~/Shared Folder Hannelore/CENTOS PROJECT HH SHARE/Making Tree in R")
cgtable <- read.csv(file = "results_alleles-brug.tsv",
                    sep = "\t", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)
metadata <- read.csv(file = "Metadata-brug.csv",
                     header = TRUE, stringsAsFactors = FALSE)


### GETTING TO KNOW THE DATA
dim(cgtable)
dim(metadata)
class(cgtable)
class(metadata)


################################################# ALTERNATIVE 7: getMinimumSpanningTree
#http://www.smartana.co.uk/mst-tutorial/MST.html

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
cgdist1 <-dist.gene(cgmatrix, method= "pairwise", pairwise.deletion=FALSE)
cgdist1
#remark first sample is not in rows but in columns (probably because otherwise empty row?)
#terug in matrix omzetten
cgdist <-as.matrix(cgdist1)
#make graph
cggraph = graph_from_adjacency_matrix(cgdist, mode="undirected", weighted=TRUE, diag=TRUE)
plot(cggraph)

cggraph$name
  
  
  #####make plot
#
cgmst2 = minimum.spanning.tree(cggraph) 
plot(cgmst2)













############################## PROBEERSELS


# Set plot layout
#layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
# Set symbols
cg.pch <- metadata$txid
levels(cg.pch) <- list("0"="158836", "1" = "1280")
cg.pch <- as.numeric(as.character(cg.pch))

# Set species colors
brewer.pal(3,"Set2")
cg.col <- metadata$txid
levels(cg.col) <- list("#66C2A5"="158836",
                         "#FC8D62"="1280")
cg.col <- as.character(cg.col)

#fancyer
plot(cgmst2, col =cg.col )



#try with ggplot
ggplot(metadata)
#TO TEST: https://corybrunson.github.io/ordr/reference/stat-biplot-spantree.html

#try with emstreeR

out <- ComputeMST(cgdist)
ggplot(data = out,
       aes(x = x, y = y,
           from = from, to = to))+
  geom_point()+
  stat_MST(colour = "red", linetype = 2)



cgnetwork= igraph_to_networkD3(cggraph)
nodes = cgnetwork$nodes #hier is staal SA32 = 1
links = cgnetwork$links #hier is staal SA32 = 0
linksValue = links$value


#artificial data
set.seed(1984)
n <- 15
c1 <- data.frame(x = rnorm(n, -0.2, sd = 0.2), y = rnorm(n, -2, sd = 0.2))
c2 <- data.frame(x = rnorm(n, -1.1, sd = 0.15), y = rnorm(n, -2, sd = 0.3))
d <- rbind(c1, c2)
d <- as.data.frame(d)

check <- ComputeMST(d)
