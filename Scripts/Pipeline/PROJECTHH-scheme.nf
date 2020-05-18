#!/usr/bin/env nextflow
 
/*
========================================================================================
                cgMLST scheme creation pipeline: PROJECTHH-scheme.nf
========================================================================================
SUMMARY
Create a cg/wgMLST scheme for your analysis based on NCBI reference assemblies based on txid OR 1 sample

INPUT: 
    * 1 fastq file 
    OR
    * txid
OUTPUT: 
    *wgMLST scheme 
    *cgMLST scheme 95%
    *cgMLST scheme 99%


MUST DEFINE
    *--fastq    : give 1 fastq file (zipped)
    *--txid     : give the txid

EXAMPLE INPUT
    *nextflow run Pipeline/PROJECTHH-scheme.nf --txid "573"                 -resume #optional
    *nextflow run Pipeline/PROJECTHH-scheme.nf --fastq "DNA123.fastq.gz"      -resume #optional

AUTHOR
Hannelore Hamerlinck <hannelore.hamerlinck@hotmail.com>

----------------------------------------------------------------------------------------
*/


// =============================  Show help message ====================================
/* 
Date creation: 15/05/2020
Show the help message when asked for help
*/

def helpMessage() {
    log.info"""
    Create a cg/wgMLST scheme for your analysis based on NCBI reference assemblies based on txid OR 1 fastq sample (zipped)

    example1:   nextflow run PROJECTHH-scheme.nf --txid "573"
    example2:   nextflow run Pipeline/PROJECTHH-scheme.nf --fastq "DNA123.fastq.gz"


    --fastq     give a fastq file that's zipped
    --help      show help message
    --krakendb  give path to kraken database
    --output    give output folder name
    --refDB     give path to place/search reference assemblies (format : "refDB-txid")
    --scheme    give path to place/search generated cg/wgMLST schemes (format: "MLST-txid")
    --txid      give a known txid number (NCBI): https://www.ncbi.nlm.nih.gov/taxonomy
    
    -resume     use to continue an analysis that was run (partly) before

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


// ===========================  Set parameters ===========================
/* 
Date creation: 15/05/2020
Set all parameters
*/
 
// Define general parameters (defaults)
params.output= "$baseDir/output-schemeMLST" //nog aanpassen met outputdir
params.refDB = "/home/hannelore/PROJECTHH/Data/refDB/"
params.scheme = "/home/hannelore/PROJECTHH/Data/cgMLSTschemes/"
params.krakendbpath = "/home/hannelore/PROJECTHH/Tools/kraken-db/16S_SILVA138_k2db/" //SILVA!
params.txid = "zero"

refDBx = params.refDB
refDBa = refDBx + "accession-lists/"




//TODO: remove final / from chosen folders

outputdir = "/home/hannelore/PROJECTHH/Scripts/testoutput"


// ===========================  Divide ===========================
/* 
Date creation: 15/05/2020
Divide the flow between txid or fastq

*/


//divide between txid given or fastq given
if (params.txid=="zero"){
    //if no txid is given
    if (params.fastq){
        //=if fastq is given

        //Set File
        Channel
             .fromPath( params.fastq)
             .ifEmpty { error "Cannot find any fastq.gz files, please try again"} //=double warning...
             .set { fastq_ch }
        
         //Make miniFastq to do Kraken
        process reduce4krakenPE {
            input:
            file(fastqFile) from fastq_ch
            output:
            file("mini-fastq") into reduced4kraken
            when:
            params.fastq
            script:
            """
            zcat ${fastqFile} | head -40000 | gzip > mini-fastq
            """
        }

        //Do kraken2 to get txid, extract the report
        process kraken2 {
            publishDir "$outputdir", mode:'copy', overwrite: true
            
            input:
            file(miniFQ) from reduced4kraken

            output:
            file("report-mini-fastq.txt") into report
            
            when:
            params.fastq

            script:
            """
            kraken2 --use-names --report report-mini-fastq.txt --threads $task.cpus \
            -db $params.krakendbpath --gzip-compressed ${miniFQ} > output.txt
            """
            //=problem local computer: working with silva we can not determinated unto species level...
        }

        //Extract txid from the report (txid with first discovered species)
        process GetTxid {
            input:
            file(rep) from report

            output:
            stdout into TXID

            when:
            params.fastq

            shell:
            '''
            column -s, -t < !{rep}| awk '$4 == "S"'| head -n 1 | cut -f5 | tr -dc '0-9'
            '''
        }

        TXID.into {txida; txidb}

        process printtxid {
            input:
            val(txid) from txidb

            output:
            stdout txidX

            exec:
            """
            Txid of given fastq file is: ${txid}
            """
        }

        txidX.view {it.trim()}

    }
    else {
        //if no fastq is given
        println "please provide either a txid OR fastq file with --txid OR --fastq"
    }
}

else {
    //if txid is given use this one
    txida = params.txid

    println """The given txid is ${txida}""".stripIndent()
}



// get accesion numbers if not yet in database
process accessions {
    publishDir "$refDBa", mode:'copy', overwrite: true
        
    input:
    val txid from txida

    output:
    file("${txid}accessions.txt") into accessionlist
    val(txid) into acctxid

    when:
    acc=file("$refDBa/${txid}accessions.txt")
    acc.exists() == false
    //enkel als de txid accessionfile niet al bestaat

    script:
    """
    esearch -db assembly -query '${txid}[txid] AND "complete genome"[filter] AND "latest refseq"[filter]' \
    | esummary | xtract -pattern DocumentSummary -element AssemblyAccession > ${txid}accessions.txt
    """
    }

process refAssembly {
    publishDir "$refDBx", mode:'copy', overwrite: true
    
    input:
    file(accesions) from accessionlist
    val(txid) from acctxid

    output:
    file("refDB-${txid}/GC*") into refDB

    //when:
    //moet niet want reeds gedifineerd door accesions!

    script:
    """
    mkdir -p refDB-${txid}
    cp ${accesions} refDB-${txid}
    cd refDB-${txid}
    bit-dl-ncbi-assemblies -w ${accesions} -f fasta -j 10
    """
}



