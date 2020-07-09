#!/usr/bin/env Rscript

# Define arguments from command line
args = commandArgs(trailingOnly=TRUE)

#set working directory
workingdir <- args[1]
setwd(workingdir)
# Define input file
cgtable <- read.csv(file = args[2],
                    sep = "\t", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)

details <- dim(cgtable)

write.csv(details, file = "TESTforR.csv")