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
*  Only for illumina paired end reads (nextera)



## AUTHOR
Hannelore Hamerlinck <hannelore.hamerlinck@hotmail.com>


----------------------------------------------------------------------------------------
*/



// TODO: 
// find better way to reach adapter file
// test for zipped or not fastqs
// nextflow.config
// what happens when introducing SE files
//  use docker container (?)-> to set + make nextflow.config file!!
'''
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''



// Show help message
def helpMessage() {
    log.info"""
    Here's some help information...
    PLEASE HELP ME
    """.stripIndent()
}


//nexflow run help PROJECTHH.nf -> you need remote resource
if (params.help){
    helpMessage()
    exit 0
}


// Define general parameters (defaults)
params.reads = "$baseDir/*_{1,2}.fastq.gz"
params.SE = false
params.PE = false
params.help = false
//params.adapters = "/opt/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
params.output= "$baseDir/output"

println """\
         PROJECTHH   chosen parameters:   
         ============================
         reads        : ${params.reads}
         add...
         """
         .stripIndent()


// Define trimming parametesr
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false
// params.phix_reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Enterobacteria_phage_phiX174_sensu_lato/all_assembly_versions/GCA_002596845.1_ASM259684v1/GCA_002596845.1_ASM259684v1_genomic.fna.gz"
params.phix_reference = "$baseDir/assets/data/GCA_002596845.1_ASM259684v1_genomic.fna.gz"




// Define all folders within output
rawfastqcdir = params.output + "/rawFastqc/"
rawmultiqcdir = params.output + "/rawMultiqc/"
trimmingdir = params.output + "/Trimmed/"



// Transform reads to files and form pairs
Channel 
    .fromFilePairs( params.reads ) //form pairs
    .ifEmpty { error "No paired files for: ${params.reads}"  } //check if empty
    .set { read_pairs_ch } //make a set


// Quality control of raw data with FASTQC & MULTIQC
process rawfastqc {
    publishDir "$rawfastqcdir", mode: 'copy', overwrite: false

    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    file("fastqc_${sample_id}") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}
    fastqc --extract -t $task.cpus -q ${reads} -o fastqc_${sample_id}
    """
}

process rawmultiqc {
    publishDir "$rawmultiqcdir", mode:'copy', overwrite: false
       
    input:
    file('*') from fastqc_ch.collect() //collect nodig om van elk paar mee te nemen en niet enkel sample1 read1&2
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 

// Trimming of data

process trimmedPE {
    publishDir "$trimmingdir", mode: 'copy', overwrite: false

    input:
    set sample_id, file(forward), file(reverse) from read_pairs_ch

    output:
    file("Trimmed_${sample_id}") into trimmed_ch

    when:
    params.PE=true

    script:
    """
    mkdir Trimmed_${sample_id}
    fastqc --extract -t $task.cpus -q ${reads} -o Trimmed_${sample_id}
    java -jar /opt/Trimmomatic-0.36/trimmomatic-0.39.jar PE -threads $task.cpus -phred33
    """
}



'''
process trimmedSE {
    publishDir "$trimmingdir", mode: 'copy', overwrite: false

    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    file("Trimmed_${sample_id}") into trimmed_ch

    when:
    params.SE=true

    script:
    """
    mkdir Trimmed_${sample_id}
    fastqc --extract -t $task.cpus -q ${reads} -o Trimmed_${sample_id}
    """
}
'''