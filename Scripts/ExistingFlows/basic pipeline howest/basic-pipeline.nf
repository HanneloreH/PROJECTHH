#!/usr/bin/env nextflow
/******************************************************************************/ 
/* To get all SRR numbers for a SRA project e.g. SRP098696 use:
   wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP098696' | grep SRR | cut -f1 -d","
/******************************************************************************/ 
/* To download test samples use:
   fastq-dump --gzip -X 100000 -O /home/phulpiau/Documents/sra/ SRR5222977 SRR5222978 SRR5222979 SRR5222980 SRR5222981 SRR5222982 SRR5222983 SRR5222984 
/******************************************************************************/ 
/* Finally run the nexflow pipeline script
   nextflow run basic-pipeline.nf --mainfolder /home/phulpiau/Documents/ --threads 4
/******************************************************************************/
/* Set input parameters
/******************************************************************************/
params.mainfolder = '/home/phulpiau/Documents/'
params.threads = '2'


/******************************************************************************/
/* Set output folders 
/******************************************************************************/
fastq_files = params.mainfolder + 'sra/*.fastq.gz'
fastqc_folder = params.mainfolder + 'fastqc/'
multiqc_folder = params.mainfolder + 'multiqc/'
// tophat2_folder = params.mainfolder + 'tophat2/'
// samtools_folder = params.mainfolder + 'samtools/'


/******************************************************************************/
/* Show info 
/******************************************************************************/
log.info """\

R N A S E Q    P I P E L I N E    
==============================
main folder  : ${params.mainfolder}
threads      : ${params.threads}
samples      : ${fastq_files}
fastqc       : ${fastqc_folder}
multiqc      : ${multiqc_folder}

"""


/******************************************************************************/
/* Start 1st process: create output folders */
/******************************************************************************/
process createFolders {

    script:
    """
    mkdir -p ${fastqc_folder}
    mkdir -p ${multiqc_folder}
    """

}


/******************************************************************************/
/* Get fastq.gz files in "fastq_ch" Channel
/* fastq_ch = Channel.fromPath( fastq_files )
/******************************************************************************/
Channel
    .fromPath( fastq_files )
    .ifEmpty { error "Cannot find any fastq.gz files: $fastq_files" }
    .set { fastq_ch }


/******************************************************************************/
/* 2nd process: run FastQC on all fastq.gz files
/******************************************************************************/
process runFastQC {
    
    publishDir fastqc_folder, mode:'copy'

    input:
    file fastq from fastq_ch

    output:
    file("fastqc_${fastq}") into fastqc_ch

    script:
    """
    mkdir fastqc_${fastq}
    fastqc -t ${params.threads} -o fastqc_${fastq} ${fastq}
    """

}


/******************************************************************************/
/* 3rd process: run MultiQC on FastQC results 
/******************************************************************************/
process runMultiQC {

    publishDir multiqc_folder, mode:'copy'

    input:
    file('fastqc*/*') from fastqc_ch

    output:
    file('multiqc_report.html')

    script:
    """
    multiqc ${fastqc_folder}
    """

}


/******************************************************************************/
/* Show result "Done" or "Oops" on completion
/******************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops something went wrong!" )
}
