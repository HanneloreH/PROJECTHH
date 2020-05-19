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



// =======================================  Show help message =============================================
/* 
- Goal: Show the help message when asked for help, give all possible parametesr
- Date creation: 15/05/2020
- Last adjusted: 19/05/2020
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
    --training  give path to training file
    --txid      give a known txid number (NCBI): https://www.ncbi.nlm.nih.gov/taxonomy
    
    -resume     use to continue an analysis that was run (partly) before

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}



// ============================================= Set Parameters ========================================
/* 
- Goal: define all default parameters and folders (can by adjusted by user in command prompt)
- Date creation: 15/05/2020
- Last adjusted: 19/05/2020
*/
 
// Define general parameters (defaults)
params.output= "$baseDir/output-schemeMLST" //nog aanpassen met outputdir
params.refDB = "/home/hannelore/PROJECTHH/Data/refDB"
params.scheme = "/home/hannelore/PROJECTHH/Data/cgMLSTschemes"
params.krakendbpath = "/home/hannelore/PROJECTHH/Tools/kraken-db/16S_SILVA138_k2db" //SILVA!
params.training= "/home/hannelore/PROJECTHH/Data/TrainingFiles"
params.txid = "zero"

// Define folders
outputdir = "/home/hannelore/PROJECTHH/Scripts/testoutput"
refDBx = params.refDB
refDBa = refDBx + "/accession-lists"
trainingx = params.training
cgMLSTx = params.scheme

//TODO: remove final / from chosen folders



// ===============================================  Get txid ===========================================
/* 
- Goal: Divide the flow based on the user choice for --txid or --fastq to get the correct txid 
  When a fastq is given kraken2 analysis is performed on a fraction of the reads (10000 = mini-fastq) to determine the txid.
- Remark: Created TXID for use as txid value, create txida for further analysis
- Date creation: 15/05/2020 (* processes were already written in april'20)
- Last adjusted: 19/05/2020
- Faced problems:
    * error message if params.txid is not given -> solved by making default value "zero" and 
    * double warning if the fastq file is empty (with 1. if-statement and 2. ifEmtpy): kept this becuase
      allows else statement if no txid and no fastq is given
    * installation of kraken and database required more memory, a new Virtual machine (CENTOS)
    * kraken2: working on local computer/VB: the sytem blocks when using the kraken database 
      FIX: use smaller silvaDB instead
            SUB-problem: working with silva we can not determinate upto species level 
            SUB-fix: adjust manually in report file for tests
    * extracting txid from kraken2 is difficult because of the dollar sign that is needed to indicate the column with awk
      and because dollar is used also by Groovy for variable interpolation
      FIX: use shell-block (using shell +  ''' in stead of """ + adding !) 
    * It was difficult to extract the txid as a value
      initial code: column -s, -t < !{report}| awk '$4 == "S"'| head -n 1 | cut -f5
      FIX: using stout as output
            SUB-problem: stdout trails "newline" to the result -> can not make a new statement
            SUB-fixes: tried several fixes (save in txt file, change into variable, work with python-language...)
            this did'nt work, eventually I changed the original statement to extract only numbers: "tr -dc '0-9'"
*/


if (params.txid=="zero"){
    //=check if no txid is given

    if (params.fastq){
        //= check if fastq is given

        //Set fastq file
        Channel
             .fromPath( params.fastq)
             .ifEmpty { error "Cannot find any fastq.gz files, please try again"} //=double warning...
             .set { fastq_ch }
        
         //Make miniFastq (10000 reads)
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
        
        //Perform kraken2 to get txid from mini-fastq
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
        }

        //Extract txid from kraken2 report
        process GetTxid {
            input:
            file(rep) from report

            output:
            stdout into gtxid

            when:
            params.fastq

            shell:
            '''
            column -s, -t < !{rep}| awk '$4 == "S"'| head -n 1 | cut -f5 | tr -dc '0-9'
            '''
        }
        
        // split channel into several for  basic value (gtxid1) + print (gtxid2) + further use (txida)
        gtxid.into {gtxid1; gtxid2; txida}
    
        //define raw value???????????????????????????????????????????????????????????????????????????
        TXID =  gtxid1.view {it.trim()}

        //Print txid to user
        process printtxid {
            input:""
            val(txid) from gtxid2

            output:
            stdout txidP

            exec:
            """
            Txid of given fastq file is: ${txid}
            """
        }
        txidP.view {it.trim()}
        
    }
    else {
        //=if no fastq is given
        println "please provide either a txid OR fastq file with --txid OR --fastq"
    }
}

else {
    //=if txid is given use this one, print txid for user
    TXID = params.txid
    println """The given txid is ${TXID}""".stripIndent()

    //write variable into channel to use in processes
    txida = Channel.value(params.txid)       
}



// =================================== Get Reference Assemblies ====================================
/*
- Goal: Get reference assemblies from NCBI via txid, 
  all proceses in this section wil run if there is no txid sepcfic accession.txt file available in the DB-folder (if-clause)
- Date creation: 18/05/2020 (* processes were already written in april'20)
- Last adjusted: 19/05/2020
- Faced problems:
        * wanted to write accessions + fastas directly to database folder, is not possible, 
          fix: first needs to save it in working dir (temporary)
        * how to put ref assembly in the right directory, no output parameter: 
          fix: go into the directory
        * problems with txid from kraken (channel) and txid from parameters input 
          fix: adjust and put both in channel
        * problems checking if accession file exists
          fix: use is.file() in clause (not assert or do not make use of channel...)
        * problems implementing if-clause for existance of accession.txt file specific for txid
          fix: ??
*/

//check if the accession exists
accFile = file("${refDBa}/${TXID}accessions.txt")
if (accFile.isFile() ) { 

    println "The accession file is available continue to scheme generation for txid ${TXID}"
    txidF1=txida
}

else {
    println "no accession file available, start acquiring reference assemblies"

    // get accession-list file from NCBI with esearch
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
        //only if txid accessionfile does not exist already, 
        //remark: double check with if-clause

        script:
        """
        esearch -db assembly -query '${txid}[txid] AND "complete genome"[filter] AND "latest refseq"[filter]' \
        | esummary | xtract -pattern DocumentSummary -element AssemblyAccession > ${txid}accessions.txt
        """
    }

    // get the fasta files matching the accession numbers and save in refDB
    process refAssembly {
        publishDir "$refDBx", mode:'copy', overwrite: true
        
        input:
        file(accesions) from accessionlist
        val(txid) from acctxid

        output:
        file("refDB-${txid}/GC*")
        val(txid) into acctxid2

        //when:
        //when is not necessary because it has already been defined by accesions!

        script:
        """
        mkdir -p refDB-${txid}
        cp ${accesions} refDB-${txid}
        cd refDB-${txid}
        bit-dl-ncbi-assemblies -w ${accesions} -f fasta -j 10
        """
    }

    // get the fasta files matching the accession numbers and save in refDB
    process unzip {
        publishDir "$refDBx/refDB-${txid}/", mode:'copy', overwrite: true
        
        input:
        val(txid) from acctxid2

        output:
        file("unzipped/GC*")
        val(txid) into acctxid3

        //when:
        //when is not necessary because it has already been defined by accesions!

        script:
        """
        mkdir unzipped
        cp $refDBx/refDB-${txid}/*.gz unzipped/
        cd unzipped
        gunzip *.gz
        """
    }
    //TO DO: look into possiblity to NOT save unzipped files (to spare memory) only in temporary files!

    // make multifasta from ref assemblies
    process multi {
        publishDir "$trainingx", mode:'copy', overwrite: true
        
        input:
        val(txid) from acctxid3

        output:
        file("multi-txid${txid}.fa")
        val(txid) into txidF1

        //when:
        //when is not necessary because it has already been defined by accesions!

        script:
        """
        cat $refDBx/refDB-${txid}/unzipped/GCF* > multi-txid${txid}.fa 
        """
    }
}



// =================================== creating wgMLST ====================================
/*
- Goal: Create a wgMLST scheme based on loci from reference assemblies
- Date creation: 19/05/2020
- Last adjusted: 19/05/2020
- Faced problems:
        * Problems trying to make processes run sequential but loose from each other (based on what's already available)
          FIX:treat cgMLST analysis as one entity, and stop the script if the scheme already exists
        * Problem: wgMSLT analysis and allele calling: time consuming step!
*/

//end script if analysis is available in the given database
cgFile = file("${cgMLSTx}/MLST-${TXID}/cgMLST/cgMLST.tsv")
if (cgFile.isFile() ) {
    println """The scheme for txid ${TXID} has already been created, stopping workflow""".stripIndent()
    exit 0
}
//TO FIX: programme exists but print message is not shown

//Split txid into 2 for processes: multi and following and make sure processes aren't started too early
txidF1.into {txid2;txidF2}

//create trainingfile with prodigal
process training {
    publishDir "$trainingx", mode:'copy', overwrite: true
    
    input:
    val(txid) from txid2

    output:
    file("txid${txid}.trn")
    val(txid) into mtxid1

    when:
    train=file("$trainingx/txid${txid}.trn")
    train.exists() == false
    //only if training file does not exist already

    script:
    """
    prodigal -i $trainingx/multi-txid${txid}.fa -p single -t txid${txid}.trn 
    """
}

//create wgMLST scheme: 
//REMARK: TIME CONSUMING (>> depending on number of genomes)
process wgMLST {
    publishDir "$cgMLSTx", mode:'copy', overwrite: true
    
    input:
    val(txid) from mtxid1

    output:
    file("MLST-${txid}/wgMLST/*")
    val(txid) into mtxid2

    script:
    """
    mkdir MLST-${txid}
    mkdir MLST-${txid}/wgMLST
    cd MLST-${txid}/wgMLST
    chewBBACA.py CreateSchema -i $refDBx/refDB-${txid}/unzipped/ \
    -o MLST-${txid}/wgMLST/schema-${txid} --cpu $task.cpus --ptf $trainingx/txid548.trn
    """
}


//trial with txid 817: bact. fragilis (only 14 ref assemblies)
//limit number of reference assemblies to use in cgMLST?

