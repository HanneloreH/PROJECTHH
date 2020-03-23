#!/usr/bin/env nextflow

/*
========================================================================================
                         cgMLST pipeline
========================================================================================
## SUMMARY
Do QC, trimming, assembly and cgMLST analysis on fastq data files (format *.fastq.gz).

## MUST DEFINE
*   --reads         :path to reads (if not in same folder)
*   --PE or --SE    :paired or single end data

## DEFAULTS
*  



## AUTHOR
Hannelore Hamerlinck <hannelore.hamerlinck@hotmail.com>


----------------------------------------------------------------------------------------
*/


// Show help message
def helpMessage() {
    """
    echo parameters
    """
}

if (params.help){
    helpMessage()
    exit 0
}

// Define all parameters
params.reads = "$baseDir/*_{1,2}.fq"  //input if not defined
params.SE = false
params.PE = false





// analysis will be done for: 
process gettingFiles {
    echo true
    
    input:
    file fastq from 

    script:
    """
    Analysis will be done for $fastq
    """
}
