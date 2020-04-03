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

## EXAMPLE INPUT
nextflow run Pipeline/PROJECTHH.nf --PE --reads "/home/hannelore/PROJECTHH/Data/RawData-KP/*_{1,2}.fastq.gz"
OPTIONAL 
-resume 
-with-docker

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
params.krakendbpath = "/home/hannelore/PROJECTHH/Tools/kraken-db/16S_SILVA138_k2db/"
params.kronataxonomy = "/home/hannelore/krona/taxonomy/"


// Define all folders within output
fastq_files = params.reads
rawfastqcdir = params.output + "/1-rawFastqc/"
rawmultiqcdir = params.output + "/2-rawMultiqc/"
trimmingdir = params.output + "/3-Trimmed/"
kraken2kronadir = params.output + "/4-kraken2-krona/"
megahitdir = params.output + "/5-Assembly-megahit/"


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
    set val(sample_id), file("Trimmed_${sample_id}/TRIM-*") into trimmedPE_reads
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
    file("Trimmed_${reads.simpleName}/TRIM-*") into trimmedSE_reads
    file("fastp.*")

    when:
    params.SE

    script:
    """
    mkdir Trimmed_${reads.simpleName}
    fastp -i "${reads}" -o Trimmed_${reads.simpleName}/TRIM-${reads}
    """
    //remark: fastp for SE has automatic adapter searching
}


//split data to be able to use the channel multiple times (for kraken and for assembly)
trimmedPE_reads.into {trimmedPE_kraken; trimmedPE_assembly}
trimmedSE_reads.into {trimmedSE_kraken; trimmedSE_assembly}


// ========================= Taxonomic sequence classification ======================
/* 
Date creation: 2/4/2020

Problems: 
- installation of kraken and database required more memory, a new Virtual machine (CENTOS)
was created.
- problems fixing the output of the PE trimmed reads into kraken2krona, I wanted to start from files BUT 
than the pipeline wants to start immediately which gives no results... 
    FIX: problem was because set was made during trimming PE reads: has been deleted

*/

process kraken2kronaPE {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true
       
    input:
    set sample_id2, file(trimfq) from trimmedPE_kraken

    output:
    file("kraken2krona-${sample_id2}/*")
     
    when:
    params.PE

    script:
    // remakr conda ngs must be activated! conda activate ngs
    """
    mkdir kraken2krona-${sample_id2}
    kraken2 --use-names --report report-${sample_id2}.txt --threads $task.cpus -db $params.krakendbpath \
    --paired --gzip-compressed ${trimfq} > kraken2krona-${sample_id2}/output.kraken
    
    cat kraken2krona-${sample_id2}/output.kraken | cut -f 2,3 > kraken2krona-${sample_id2}/results.krona
    
    ktImportTaxonomy kraken2krona-${sample_id2}/results.krona -o kraken2krona-${sample_id2}/krona.html \
    --taxonomy $params.kronataxonomy
    """
}


process kraken2kronaSE {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true
       
    input:
    file(trimfq) from trimmedSE_kraken

    output:
    file("kraken2krona-${trimfq.simpleName}/*")
    
    when:
    params.SE

    script:
    """
    mkdir kraken2krona-${trimfq.simpleName}
    kraken2 --use-names --report report-${trimfq.simpleName}.txt --threads $task.cpus -db $params.krakendbpath \
    --gzip-compressed ${trimfq} > kraken2krona-${trimfq.simpleName}/output.kraken
    
    cat kraken2krona-${trimfq.simpleName}/output.kraken | cut -f 2,3 > kraken2krona-${trimfq.simpleName}/results.krona
    
    ktImportTaxonomy kraken2krona-${trimfq.simpleName}/results.krona --taxonomy $params.kronataxonomy
    """
}



// ================================= ASSEMBLY  =============================
/* 
Date creation: 3/4/2020

Problems: /

Possibilities: include other assemblers and give a choice to users

*/

process megahitPE {
    publishDir "$megahitdir", mode:'copy', overwrite: true
       
    input:
    set sample_id, file(trimfq) from trimmedPE_assembly

    output:
    file("megahit-${sample_id}/*") into megahitPE
    file("megahit-${sample_id}/*") into trimmedSE_reads
     
    when:
    params.PE

    script:
    """
    megahit -1 "${trimfq[0]}" -2 "${trimfq[1]}" -o megahit-${sample_id}
    """
}

process megahitSE {
    publishDir "$megahitdir", mode:'copy', overwrite: true
       
    input:
    file(trimfq) from trimmedSE_assembly

    output:
    file("megahit-${trimfq.simpleName}/*")
     
    when:
    params.SE

    script:
    """
    megahit -r "${trimfq}" -o megahit-${trimfq.simpleName}
    """
}







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
// fix installation krona, via conda?? conda activate ngs
    // also fastp via bioconda -> problems??
// fix krona: taxnoomy db (so no 0 result :p)
// make a list of short names (for long fastq names)
// check if fastp.html gives view of all samples


'''
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''
