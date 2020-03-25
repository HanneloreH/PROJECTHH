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

'''
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
'''
// TODO: find better way to reach adapter file


// Define all parameters (defaults)
params.reads = "$baseDir/*_{1,2}.fq"
params.SE = false
params.PE = false
params.adapters = "/opt/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
params.output= "$baseDir/output"


// Define files/paths/...
reads_files = file(params.reads)


// analysis will be done for: 
process Pfastqc {
     //publishDir params.output, mode: 'copy', overwrite: false
    
    input:
    file fastq from reads_files

    output:
    file("output") into fastqc_ch

    script:
    """
    mkdir output
    fastqc --extract -t $task.cpus -q ${fastq} -o output/
    """
}









''''
// default use docker container -> to set + make nextflow.config file!!
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''


'''

//if PE is true form pairs
read_pairs_ch = Channel .fromFilePairs(params.reads)




// analysis will be done for: 
process gettingFiles {
     publishDir "$baseDir/chunks", mode: 'copy', overwrite: false
    
    input:
    file fastq from 

    script:
    fastqc --extract -t $taks.cpus -o 
}

To acces number of cpu's
(--threads) $task.cpus

'''

