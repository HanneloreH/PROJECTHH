#!/usr/bin/env Rscript


### INSTALL & LOAD LIBRARIES
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ape, emstreeR, ggplot2, ggrepel, reshape2)
#library(ape)
#library(emstreeR)
#library(ggplot2)
#library(ggrepel)
#library(reshape2)



### DEFINE DATA
# Define arguments from command line
args = commandArgs(trailingOnly=TRUE)
#set working directory (where files will be saved): DOES NOT WORK (no sudo)
#workingdir <- args[1]
#setwd(workingdir)
# Define input file
cgtable <- read.csv(file = args[1],
                    sep = "\t", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)
metadata <- read.csv(file = args[2],
                    sep = ",", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)
resultdir <- args[3]





### CLEAN UP DATA
#change rownames (from numbers to actual sample names)
row.names(cgtable)<-cgtable$FILE
#remove first column (with samples names)
cgtable<-cgtable[,-1]




### PREPPING FOR PLOT
## 1)Distance matrix (pairwise, alleles)
#library(ape)
cgmatrix <- as.matrix(cgtable)
cgdist <-dist.gene(cgmatrix, method= "pairwise", pairwise.deletion=FALSE)
cgdistMAT <- as.matrix(cgdist)
#save the file at resultdir
distname = paste(resultdir, 'DistanceMatrix.csv',sep= '')
write.csv(cgdistMAT, file = distname)

## 2)minimum spanning tree
#library(emstreeR)
#create mst from distance matrix
cgmst<-ComputeMST(as.matrix(cgdist))

## 3) Principal component analysis
#princ. comp analysis
PrinC <- prcomp(cgdistMAT)




### PLOTTING
#label names !!! make sure samples in metadata are in same order!!!
Names = metadata[,1]

#ggplot for every metadata column:
for (i in 1:ncol(metadata)){

    #set colours
  cgmst$coloring <-as.factor(metadata[,i])
  title <-colnames(metadata)[i]
  
  
  #plot
  ggplot(data = cgmst, 
         #plot based on principle coordinate points
         aes(x = PrinC$x[, 1], y = PrinC$x[, 2], 
             from = from, to = to, colour=coloring))+ 
    labs(color=title)+
    #make sure data points are not overlapping with position_dodge2
    geom_point(size=3,position= position_dodge2(width=100))+
    # geom_text(aes(label=Names),hjust=0,vjust=1)+
    stat_MST(colour = "black", linetype = 1)+ 
    geom_label_repel(aes(label = Names),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50')+
    labs(title="Minimum Spanning Tree cgMLST data", subtitle="OUTB8", 
         x="absolute number of allele differences", y="absolute number of allele differences")
  
  #define file name
  ctitle <-colnames(metadata)[i]
  filename <- paste0(resultdir,"/MST-plot-color-", ctitle,".jpeg")
  
  #save file
  ggsave(filename)
  
}
