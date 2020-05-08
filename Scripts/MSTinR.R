#Making a minimal spanning tree in R


### LOAD REQUIRED PACKAGES/LIBRARIES
library("igraph") #for MST



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

## MAKE DISTANCE MATRIX
cgDist<-dist(cgtable, method="euclidean", diag = FALSE, upper = FALSE, p = 2)
dim(cgDist)

#creating an igraph graph from a data frame
cgGraph <- graph_from_data_frame(cgtable, directed=TRUE, vertices=NULL)
class (cgGraph)





mst(cgGraph, weights=NULL, algorithm=NULL)





################################################ ALTERNATIVE

#INSTALL

library(curl)
library(quantmod)
library(devtools)

#WERKT NIET...
devtools::install_github('systematicinvestor/SIT.date')
curl_download('https://github.com/systematicinvestor/SIT/raw/master/SIT.tar.gz', 'sit',mode = 'wb',quiet=T)
install.packages('sit', repos = NULL, type='source')
#installatie uit tar.gz file BUT then "SIT.dat is not avilable'
library(sit)

