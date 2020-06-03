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
Last adjusted: 8/04/2020 (folders)

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
params.krakendbpathsilva = "/home/hannelore/PROJECTHH/Tools/kraken-db/16S_SILVA138_k2db/"
params.krakendbpath = "/home/hannelore/PROJECTHH/Tools/kraken-db/minikraken_8GB_202003/minikraken_8GB_20200312"
params.kronataxonomy = "/home/hannelore/krona/taxonomy/"
params.refDB = "/home/hannelore/PROJECTHH/Data/refDB"

// Define all folders within output
fastq_files = params.reads
qualitydir = params.output + "/1-Quality/"
trimmingdir = params.output + "/2-Trimmed&Quality/"
reduced4krakendir = params.output + "/3-kraken2-krona/reduced"
kraken2kronadir = params.output + "/3-kraken2-krona/"
megahitdir = params.output + "/4-Assembly-megahit/"
quastdir = params.output + "/4-Assembly-megahit/quast/"

// Folders associated with pipeline
refAssemblyDBdir = params.refDB
refAssemblyACCdir = params.refDB + "/accesion-lists/"
krakendbpath=params.krakendbpath

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
Last adjustment: 8/04/2020 (adjusted output)

Remarks: do fastqc, make one folder per sample (if PE: both FWD as REV)

Problems during creation: only one file wad done due to wrong input during run of the command,
"" were not added

Todo: 
*/

// Quality control of raw data with FASTQC
process rawfastqc {
    publishDir "$qualitydir", mode: 'copy', overwrite: true

    input:
    file fastq from fastq_ch

    output:
    file("fastqc_${fastq}/*.zip") into fastqc_ch
    file("fastqc_${fastq}/*.html")

    script:
    """
    mkdir fastqc_${fastq}
    fastqc --extract -t $task.cpus -q ${fastq} -o fastqc_${fastq}
    """
}

// do MULTIQC
process rawmultiqc {
    publishDir "$qualitydir", mode:'copy', overwrite: true
       
    input:
    file('*') from fastqc_ch.collect() 
    //collect nodig om van elk paar mee te nemen en niet enkel sample1 read1&2
    
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
Last adjusted: 8/04/2020

Define with kraken what your species txid is in the sample (needed for scaffolding / cgMLST allele calling)
to lower computational needs do analysis on mini-fastq

Problems: 
- installation of kraken and database required more memory, a new Virtual machine (CENTOS)
was created.
- problems fixing the output of the PE trimmed reads into kraken2krona, I wanted to start from files BUT 
than the pipeline wants to start immediately which gives no results... 
    FIX: problem was because set was made during trimming PE reads: has been deleted

*/


//FIRST take 10000 reads to test with kraken (to limit needed time for analysis)
process reduce4krakenPE {
    publishDir "$reduced4krakendir ", mode:'copy', overwrite: true
       
    input:
    set sample_id, file(trimfq) from trimmedPE_kraken

    output:
     set val(sample_id), file("reduced-*") into reducedPE_kraken
     
    when:
    params.PE

    script:
    """
    zcat ${trimfq[0]} | head -40000 | gzip > reduced-${trimfq[0]}
    zcat ${trimfq[1]} | head -40000 | gzip > reduced-${trimfq[1]}
    """
}
process reduce4krakenSE {
    publishDir "$reduced4krakendir ", mode:'copy', overwrite: true
       
    input:
    file(trimfq) from trimmedSE_kraken

    output:
    file("reduced-*") into reducedSE_kraken
     
    when:
    params.SE

    script:
    """
    zcat ${trimfq} | head -40000 | gzip > reduced-${trimfq}
    """
}

// SECOND: do kraken and prepare for krona
process kraken2kronaPE {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true
       
    input:
    set sample_id2, file(trimfq) from reducedPE_kraken

    output:
    file("kraken2krona/output-${sample_id2}.kraken")
    file("kraken2krona/report-${sample_id2}.txt") into krakenPE4txid
    file("kraken2krona/results-${sample_id2}.krona") into krakenPE4krona
    
    when:
    params.PE

    script:
    """
    mkdir kraken2krona
    kraken2 --use-names --report kraken2krona/report-${sample_id2}.txt \
    --threads $task.cpus -db $params.krakendbpathsilva \
    --paired --gzip-compressed ${trimfq} > kraken2krona/output-${sample_id2}.kraken
    
    cat kraken2krona/report-${sample_id2}.txt | cut -f 2,3 > kraken2krona/results-${sample_id2}.krona
    """
}

//ONLY ON SERVER: (use 8Gb database)
'''
process kraken2kronaPE8Gb {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true
       
    input:
    set sample_id2, file(trimfq) from reducedPE_kraken

    output:
    file("kraken2krona/output-${sample_id2}.kraken")
    file("kraken2krona/report-${sample_id2}.txt") into krakenPE4txid
    file("kraken2krona/results-${sample_id2}.krona") into krakenPE4krona
    
    when:
    params.PE

    script:
    """
    mkdir kraken2krona
    kraken2 --use-names --report kraken2krona/report-${sample_id2}.txt \
    --threads 1 -db $krakendbpath \
    --paired --gzip-compressed ${trimfq} > kraken2krona/output-${sample_id2}.kraken
    
    cat kraken2krona/report-${sample_id2}.txt | cut -f 2,3 > kraken2krona/results-${sample_id2}.krona
    """
}
'''

process kraken2kronaSE {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true
       
    input:
    file(trimfq) from reducedSE_kraken

    output:
    file("kraken2krona/output-${trimfq.simpleName}.kraken")
    file("kraken2krona/report-${trimfq.simpleName}.txt") into krakenSE4txid
    file("kraken2krona/results-${trimfq.simpleName}.krona") into krakenSE4krona
    file("kraken2krona-/txid-${trimfq.simpleName}") into txidSE //werkt nog niet!!
    
    when:
    params.SE

    script:
    """
    mkdir kraken2krona
    kraken2 --use-names --report kraken2krona/report-${trimfq.simpleName}.txt --threads $task.cpus \
    -db $params.krakendbpathsilva --gzip-compressed ${trimfq} > kraken2krona/output${trimfq.simpleName}.kraken
    
    column -s, -t < kraken2krona/report-${trimfq.simpleName}.txt | awk '!{4} == "S"'| head -n 1 | cut -f5 \
    > kraken2krona-/txid-${trimfq.simpleName}

    cat kraken2krona/report-${trimfq.simpleName}.txt | cut -f 2,3 > kraken2krona/results-${trimfq.simpleName}.krona
    """
}

// THIRD: krona
// remark conda ngs must be activated for krona! "conda activate ngs"
// update taxonomy
process krona_db {
    output:
    file("taxonomy/taxonomy.tab") into krona_taxonomy

    script:
    """
    ktUpdateTaxonomy.sh taxonomy
    """
}

// create mutual channel for PE and SE
report4krona = Channel.create()
if (params.PE) {
report4krona = krakenPE4krona
}
else {
report4krona = krakenSE4krona
}

process krona {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true
    input:
    file("taxonomy/taxonomy.tab") from krona_taxonomy
    file (report) from report4krona

    output:
    file("*.html")

    script:
    """
    ktImportTaxonomy ${report} -tax taxonomy -o KRONA-${report.simpleName}.html

    """
}


// ================================= REFERENCE ASSEMBLIES  =============================
/* 
Date creation: 9/4/2020
Last modification: 15/4/2020 (introducing python to get txid)

Get reference assemblies from NCBI for the given txid, necessary for scaffolding AND cgMLST scheme creation

Problems/fixes:
- takes a long time to download, would be good to have a database (txid)
- problem extracting txid from kraken2 because of dollar sign needed to indicate column with awk: 
    * fix: shell-block(using shell +  ''' in stead of """ + adding !) to fix the problem of $4 : 
            because dollar is used also by Groovy for variable interpolation
- problem 2 how to get txid in a value: 
    FIX
        * BASH first (shell-block): column -s, -t < !{report}| awk '$4 == "S"'| head -n 1 | cut -f5
        * BASH better: awk -F, '$4 == "S" { printf("%s", $5); exit }' "!{report}"
        * HOW OUTPUT?
            - Fixed by using stdout as output
            - New problem: stdout trails "newline" to the result -> can not make a new statement!
                *Fix2.1: make .txt: can't find way to read it without stout: not found: FAIL
                *Fix2.2: change into variable: not found: FAIL
                *Fix2.3: change language: python!
                    PROBLEM: variable is also not working
                    Save script seperately as pythong module: FAIL, module will not load
                        see ReadTxid.py (in folder Tests)
                *Fix2.4: write groovy script to get the number: FAIL
                *Fix2.5: back to bash: add tr -dc '0-9' to extract only numbers
- problem3: wanted to write accessions + fastas directly to database, is not possible, first needs to go to working dir!
- Problem4: how to put ref assembly in the right directory, no output parameter: FIX: go into the directory

//remark with silva database highest level is Genus (G) in stead of (S) (more correct)!!!!

*/

// create mutual report channel for PE and SE
report4txid = Channel.create()
if (params.PE) {
report4txid = krakenPE4txid
}
else {
report4txid = krakenSE4txid
}

// extract the txid with bash: too many problems -> go 
process txid {
    publishDir "$kraken2kronadir", mode:'copy', overwrite: true

    input:
    file(report) from report4txid

    output:
    stdout into txid4assembly

    shell:
    '''
    column -s, -t < !{report}| awk '$4 == "S"'| head -n 1 | cut -f5 | tr -dc '0-9'
    '''
    //IF SILVA DB: "G" (OPM: manually adjusted the trial files with silva database to S: 573)
    //IF 8GB DB :"S"
}

txid4assembly.into {txid4assembly_print; txid4assembly_accesion; txid4assembly_exists}

//print txid - OPTIONAL
process printtxid {

    input:
    val(txid) from txid4assembly_print

    output:
    stdout result

    exec:
    """
    echo txid of reference assembly is: ${txid}
    """
    //use get.text() or print?
}
result.view {it.trim()}

// get accesion numbers if not yet in database
process accessions {
    publishDir "$refAssemblyACCdir", mode:'copy', overwrite: true
    
    input:
    val txid from txid4assembly_accesion

    output:
    file("${txid}accessions.txt") into accessionlist
    val(txid) into accesontxid

    when:
    acc=file("$refAssemblyACCdir/${txid}accessions.txt")
    acc.exists() == false
    //enkel als de txid accessionfile niet al bestaat

    script:
    """
    esearch -db assembly -query '${txid}[txid] AND "complete genome"[filter] AND "latest refseq"[filter]' \
    | esummary | xtract -pattern DocumentSummary -element AssemblyAccession > ${txid}accessions.txt
    """
}

process refAssembly {
    publishDir "$refAssemblyDBdir", mode:'copy', overwrite: true
    
    input:
    file(accesions) from accessionlist
    val(txid) from accesontxid

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

'''
//NIET NODIG: ALLEEN PATH
process fetchRefAssembly {
    
    input:
    val(txid) from txid4assembly_exists

    output:


    when:
    //enkel als de ref database WEL al bestaat
    acc=file("$refAssemblyDBdir/refDB-${txid}/GC*")
    acc.exists() == true

    script:
    """
    cd $refAssemblyDBdir/refDB-${txid}

    """
}
'''



// ================================= ASSEMBLY  =============================
/* 
Date creation: 3/4/2020

Make de novo assembly with megahit and QA of the assembly with quast

Problems: 
- problems with matplotlib (quast)
- output assembly = final.contigs.fa, but then there's no samplename anymore!

Possibilities: include other assemblers and give a choice to users

*/

//import matpotlib
process matplotlib {
    script:
    '''
    #!/usr/bin/env python3
    import matplotlib
    '''
}

//perform assembly with megahit (PE and SE)
process megahitPE {
    publishDir "$megahitdir", mode:'copy', overwrite: true
       
    input:
    set sample_id, file(trimfq) from trimmedPE_assembly

    output:
    file("megahit/*.fa") into megahitPE    
    
    when:
    params.PE

    script:
    """
    megahit -1 "${trimfq[0]}" -2 "${trimfq[1]}" -o megahit \
    --out-prefix ${sample_id} -t $task.cpus
    """
}
process megahitSE {
    publishDir "$megahitdir", mode:'copy', overwrite: true
       
    input:
    set val(sample_id), file("megahit/*.fa") from trimmedSE_assembly

    output:
    file("megahit/*.fa") into megahitSE 
     
    when:
    params.SE

    script:
    """
    megahit -r "${trimfq}" -o megahit \
    --out-prefix ${trimfq.simpleName} -t $task.cpus
    """
}

// split the file and stop division in SE and PE for chewbbaca
if (params.PE){
    megahitPE.into {megahit_quast; megahit_chew}
}
else {
    megahitSE.into {megahit_quast; megahit_chew}
}


// Control of assembly with Quast
process quast{
    publishDir "$quastdir", mode:'copy', overwrite: true
       
    input:
    file(assembly) from megahit_quast

    output:
    file("quast-${assembly.simpleName}/*")

    script:
    """
    mkdir quast-${assembly.simpleName}/
    quast.py ${assembly} -o quast-${assembly.simpleName}/
    """
}


'''
// ================================= Scaffolding =============================
/* 
Date creation: 3/4/2020

Problems: 


Possibilities:
*/


val txid from txid4assembly_exists 
    //because assembly takes long time it's not a problem for samples which do not yet have a
    //ref before starting the protoco




// ================================= cgMLST analysis =============================
/* 
Date creation: 3/4/2020

Problems: 


Possibilities:

*/


'''





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
// database with already downloaded assemblies per txid


'''
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''
