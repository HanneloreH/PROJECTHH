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
*  Only for illumina  paired end reads (nextera)

## AUTHOR
Hannelore Hamerlinck <hannelore.hamerlinck@hotmail.com>

----------------------------------------------------------------------------------------
*/



// =============================  Show help message ====================================
/* 
Date creation: 25/03/2020

Remarks: first the helpmessage is defined in a function. Standard params.help=false (also 
see nextflow.contig) but when used in command this is set to true and the function is
activated with an if-clause.

Problems during creation: in the beginning this set-up did not print the help-message, 
this was overcome by adding "log.info"

Todo: add help-message itself :)
*/

def helpMessage() {
    log.info"""
    Here's some help information...
    PLEASE HELP ME
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


// ===========================  Set parameters ===========================
/* 
Date creation: 25/03/2020

Remarks: all parameters (files, settings, directories) are put together to find/adjust
them easily

Todo: add extra when necessary
*/
 
// Define general parameters (defaults)
params.reads = "$baseDir/*_{1,2}.fastq.gz"
params.SE = false
params.PE = false
params.help = false
params.adapters = "/opt/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
params.output= "$baseDir/output"
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false


// Define all folders within output
fastq_files = params.reads
rawfastqcdir = params.output + "/rawFastqc/"
rawmultiqcdir = params.output + "/rawMultiqc/"
trimmingdir = params.output + "/Trimmed/"


// Transform reads to files and form pairs and set channels
Channel
    .fromPath( params.reads)
    .ifEmpty { error "Cannot find any fastq.gz files: $fastq_files" }
    .set { fastq_ch }

//for PE files
Channel 
    .fromFilePairs( params.reads ) //form pairs
    .ifEmpty { error "No paired files for: ${params.reads}"  } //check if empty
    .set { read_pairsPE_ch } //make a set

//For SE files
Channel
    .fromPath( params.reads ) //get data from path
    .ifEmpty { error "Cannot find any files: ${params.reads}" } //check if empty
    .set { read_pairsSE_ch } //make a set
'''
if (params.SE){
}
else {
    read_pairsSE_ch=false
}
'''

// Print the relevant parameters set
println """\
         PROJECTHH   chosen parameters:   
         ============================
         reads        : ${params.reads}
         fastq_files      : ${fastq_files}
         channels:
         - fastq_ch : ${fastq_ch}
         """
         .stripIndent()




// ===========================  Quality control raw data ===========================
/* 
Date creation: 25/03/2020

Remarks: do fastqc, make one folder per sample (if PE: both FWD as REV)

Problems during creation: only one file wad done due to wrong input during run of the command,
"" were not added

Todo: 
*/

// Quality control of raw data with FASTQC
process rawfastqc {
    publishDir "$rawfastqcdir", mode: 'copy', overwrite: true

    input:
    file fastq from fastq_ch

    output:
    file("fastqc_${fastq}") into fastqc_ch

    script:
    """
    mkdir fastqc_${fastq}
    fastqc --extract -t $task.cpus -q ${fastq} -o fastqc_${fastq}
    """
}

// do MULTIQC
process rawmultiqc {
    publishDir "$rawmultiqcdir", mode:'copy', overwrite: true
       
    input:
    file('*') from fastqc_ch.collect() //collect nodig om van elk paar mee te nemen en niet enkel sample1 read1&2
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 


// ===========================  Data Trimming ===========================
/* 
Date creation: 27/03/2020

Remarks: Chosen fastp for trimming + adaptor clipping because easier in use than trimmomatic
+ faster + direct quality parametesr

Todo: add extra parametesr to fastp
    -g for polyG tail trimming, -x for polyX tail trimming
    -p for overrepresented sequence analyses, takes time!
    --detect_adapter_for_pe
*/


// Trimming of PE data
process trimmedPE {
    publishDir "$trimmingdir", mode: 'copy', overwrite: true

    input:
    set sample_id, file(reads) from read_pairsPE_ch

    output:
    set val(sample_id), file("Trimmed_${sample_id}/TRIM-*.fastq.gz") into trimmedPE_reads
    file("fastp.*")

    when:
    params.PE

    script:
    """
    mkdir Trimmed_${sample_id}
    fastp --detect_adapter_for_pe -i "${reads[0]}" -I "${reads[1]}" \
    -o Trimmed_${sample_id}/TRIM-${reads[0]} -O Trimmed_${sample_id}/TRIM-${reads[1]}
    """
}


// Trimming of SE data
process trimmedSE {
    publishDir "$trimmingdir", mode: 'copy', overwrite: false

    input:
    file(reads) from read_pairsSE_ch

    output:
   file("Trimmed_${reads}/TRIM-*.fastq.gz") into trimmedSE_reads
    file("fastp.*")

    when:
    params.SE

    script:
    """
    mkdir Trimmed_${reads}
    fastp --detect_adapter_for_pe -i "${reads[0]}" -o Trimmed_${reads}/TRIM-${reads[0]}
    """
}


// ========================= Taxonomic sequence classification ======================
/* 
Date creation: 27/03/2020

use: kraken and krona
*/


process rawmultiqc {
    publishDir "$rawmultiqcdir", mode:'copy', overwrite: true
       
    input:
    file('*') from fastqc_ch.collect() //collect nodig om van elk paar mee te nemen en niet enkel sample1 read1&2
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """




// =============================  Finished message ====================================
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops something went wrong!" )
}




// TODO: 
// find better way to reach adapter file
// test for zipped or not fastqs
// nextflow.config
// fastqc en multiqc after trimming
//  use docker container (?)-> to set + make nextflow.config file!!


'''
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''
