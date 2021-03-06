#!/usr/bin/env nextflow
 
/*
========================================================================================
                         cgMLST analysis on fastq data
========================================================================================
SUMMARY
Do WGS bacterial analysis on fastq files (format *.fastq.gz) based on a known scheme for cg/wgMLST (includes trimming and assembly)

INPUT (must define):
    * --PE or --SE  :paired or single end data
    * --reads       :path to reads (if not in same folder)
    * --scheme      :path to scheme (folder with fasta files defined by scheme)

OUTPUT: 
    * QA results before and after trimming in folder "Quality" 
    * Assembly for every sample and assembly quality parameters in folder "Assemblies-*"
    * Minimum Spanning Tree and conclusion in folder "Results"

EXAMPLE INPUT:
    * nextflow run Pipeline/PROJECTHH-analysis2.nf --PE --reads "/home/hannelore/PROJECTHH/Data/RawData-KP/" \
      --scheme "/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-287-all/cgMLST/scheme-287-all-cgMLST"


AUTHOR
Hannelore Hamerlinck <hannelore.hamerlinck@hotmail.com>
----------------------------------------------------------------------------------------
*/



// =============================  Show help message ====================================
/* 
- Date creation: 25/03/2020
- Last adjusted: 3/06/2020
- Goal: Give help message if asked
- Remarks: first the helpmessage is defined in a function. Standard params.help=false (also 
  see nextflow.contig) but when used in command this is set to true and the function is
  activated with an if-clause.
- Faced problems:
    * in the beginning this set-up did not print the help-message,this was overcome by adding "log.info"
*/

def helpMessage() {
    log.info"""
    Do WGS bacterial analysis on fastq files based on a known scheme for cg/wgMLST

    example:   nextflow run Pipeline/PROJECTHH.nf --PE --reads "Data/*_{1,2}.fastq.gz" --scheme "Scheme/cgMLST"

    --assem     OPTIONAL give path to established assemblies that need to be included in the analysis (format: *.fa)
    --cpu       give maximal number of CPUs (default = 1)
    --help      show help message
    --output    give path to output folder
    --PE        use for paired end data
    --reads     give path to input fastq files (format *.fastq.gz)
    --scheme    give path to wg/cg scheme folder 
    --SE        use for single end data
    --training  give path to training file 
    
    -resume     use to continue an analysis that was run (partly) before
    -with-docker ***under construction***

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}



// ===========================  Set parameters ===========================
/* 
- Date creation: 25/03/2020
- Last adjusted: 3/06/2020
- Goal: All parameters (files, settings, directories) are put together to find/adjust
  them easily
- Faced problems:
    *Problem with input paths: some have "/" at the end others did not, fixed by adding a function that trims the final /
*/
 
// Define general parameters (defaults)
params.assem = false
params.cpu= 1
params.help = false
params.output= "$baseDir/output"
params.PE = false
params.reads = "$baseDir/*{1,2}.fastq.gz"
params.scheme = "$baseDir/cgMLST"
params.SE = false
params.training = "$baseDir/*.trn"

// Define all folders 
//Remove final "/" from workingfolders and define them:
//function to trim final /
def trimFolder = { 
    it.endsWith("/") ? it[0..-2] : it
}
//input folders
fastq_f = trimFolder("$params.reads")
fastq_files = fastq_f + "/*{1,2}.fastq.gz"
inputscheme = trimFolder("$params.scheme")
trainingfile = trimFolder("$params.training")
assdir = trimFolder("$params.assem")
assemdir = assdir + "/*.fa"
//output folders
outputdir = trimFolder("$params.output")
qualitydirPRE = outputdir + "/Quality/raw"
qualitydirPOST = outputdir + "/Quality/trimmed"
megahitdir = outputdir + "/Assemblies-megahit"
mlstdir = outputdir + "/MLSTtypes"
quastdir = megahitdir + "/assembly-quality"
analysisdir = outputdir + "/Analysis"
input4dir = analysisdir + "/input4analysis"
//resultdir = outputdir + "/Results"  //defined in analysis

// Transform reads to files and form pairs and set channels
//1) Input channel for quality control
Channel
    .fromPath("${fastq_files}")
    .ifEmpty { error "Cannot find any fastq.gz files: $fastq_files" }
    .set {fastq_ch}
//2) Input channel for further analysis (split between PE and SE reads with if-clause)
if (params.SE){
    Channel
    .fromPath("${fastq_files}") //get data from path
    .ifEmpty { error "Cannot find any files: $fastq_files" } //check if empty
    .set {read_pairs_ch} //make a set
}
else {
    Channel 
    .fromFilePairs("${fastq_files}") //form pairs
    .ifEmpty { error "No paired files for: $fastq_files"  } //check if empty
    .set {read_pairs_ch} //make a set
}

// Print the parameters set
println """\
        Analysis starting using following parameters:   
        ==============================================================
        * Path to extra assemblies        : ${assemdir}
        * Number of CPUs                  : ${params.cpu}
        * Output-folder                   : ${outputdir}
        * Folder with input reads         : ${fastq_f}
        * Folder with chosen MLST scheme  : ${inputscheme}
        * Training file to use            : ${trainingfile} 
        * PE reads                        : ${params.PE}
        * SE reads                        : ${params.SE}
         """
         .stripIndent()


//SEQUENTIAL
//Make sure everything happens sequential
SEQ1 = Channel.value( 'x' )
SEQ1.into {seq1a; seq1b}


// ===========================  Quality control raw data ===========================
/* 
- Date creation: 25/03/2020
- Last adjustment: 3/6/2020
- Goal: Do multiqc on raw data
- Remarks: one file is made per sample (if PE: both FWD as REV)
- Faced problems:
    * only one file was done (of batch) due to wrong input at command line, this was fixed with ""
    * Not all output should be copied with publishdir, fixed with "pattern"
*/

// Quality control of raw data with FASTQC
process rawfastqc {
    publishDir "$qualitydirPRE", mode: 'copy', overwrite: true, pattern: "fastqc_${fastq}/*.html"

    input:
    file fastq from fastq_ch

    output:
    file("fastqc_${fastq}/*.zip") into fastqc_ch
    file("fastqc_${fastq}/*.html")

    script:
    """
    mkdir fastqc_${fastq}
    fastqc --extract -t $params.cpu -q ${fastq} -o fastqc_${fastq}
    """
}

// do MULTIQC
process rawmultiqc {
    publishDir "$qualitydirPRE", mode:'copy', overwrite: true
       
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
- Date creation: 27/03/2020
- Last adjustment: 3/6/2020
- Goal: Trim the raw files
- Remarks: I chose "fastp" for trimming + adaptor clipping because easier in use than trimmomatic
  + faster + direct quality parametesr
- Faced problems:
    * Make command that takes into account adding a certain parameter or not
      FIX: ....
    * Originally I made two different processes (one for PE and one for SE data), use was defined by when-parameter
      It was pointed to me that this may be combined into one process to remove redundant code
      I tried this using an if-clause at several positions in the process "trimming", 
      however this gave errors which I was not able to fix, that's why I remained keeping 2 split 
      processes for as long as necessary.

Todo: add extra parametesr to fastp
    -g for polyG tail trimming, -x for polyX tail trimming
    -p for overrepresented sequence analyses, takes time!
    --detect_adapter_for_pe
*/

 //Trimming
 //split file
 read_pairs_ch.into {read_pairsSE_ch; read_pairsPE_ch}

 if (params.SE){
    // Trimming of SE data
    process trimmedSE {
        publishDir "$qualitydirPOST", mode: 'copy', overwrite: false,  pattern: "fastp.html"

        input:
        file(reads) from read_pairsSE_ch
        val(seq) from seq1a

        output:
        file("Trimmed_${reads.simpleName}/TRIM-*") into trimmedSE_reads
        file("Trimmed_${reads.simpleName}/TRIM-*") into trimmedSE_quality
        file("fastp.html")
        val(seq) into seq2a


        when:
        params.SE

        script:
        """
        mkdir Trimmed_${reads.simpleName}
        fastp -i "${reads}" -o Trimmed_${reads.simpleName}/TRIM-${reads}-${seq}
        """
        //remark: fastp for SE has automatic adapter searching
    }
 }
 else{
    //Trimming of PE data
    process trimmedPE {
        publishDir "$qualitydirPOST", mode: 'copy', overwrite: true, pattern: "fastp.html"

        input:
        set sample_id, file(reads) from read_pairsPE_ch
        val(seq) from seq1b

        output:
        set val(sample_id), file("Trimmed_${sample_id}/TRIM-*fastq.gz") into trimmedPE_reads
        file("Trimmed_${sample_id}/TRIM-*fastq.gz") into trimmedPE_quality
        file("fastp.html")
        val(seq) into seq2b

        when:
        params.PE

        script:
        """
        mkdir Trimmed_${sample_id}
        fastp --detect_adapter_for_pe -i "${reads[0]}" -I "${reads[1]}" \
        -o Trimmed_${sample_id}/TRIM-${reads[0]}-${seq} -O Trimmed_${sample_id}/TRIM-${reads[1]}-${seq}
        """
    }
 }
// Quality control on PE trimmed data
 if (params.SE){
    //Quality control on SE trimmed data
    process trimfastqcSE {

        input:
        file (fastq) from trimmedSE_quality
        val(seq) from seq2a

        output:
        file("*.zip") into fastqc_SE_trim_ch
        file("*.html")
        val(seq) into seq3a

        script:
        """
        fastqc --extract -t $params.cpu -q ${fastq}
        """
    }
 }
else{
    process trimfastqcPE {

        input:
        file (fastq) from trimmedPE_quality
         val(seq) from seq2b

        output:
        file("*.zip") into fastqc_PE_trim_ch
        file("*.html")
        val(seq) into seq3b

        script:
        """
        fastqc --extract -t $params.cpu -q ${fastq[0]}
        fastqc --extract -t $params.cpu -q ${fastq[1]}
        """
    }
}
//merge into one channel
if (params.PE){
    fastqc_trim_ch=fastqc_PE_trim_ch
}
else {
    fastqc_trim_ch=fastqc_SE_trim_ch
}
// do MULTIQC after trimming
process trimmultiqc {
    publishDir "$qualitydirPOST", mode:'copy', overwrite: true
       
    input:
    file('*') from fastqc_trim_ch.collect() 
    //collect nodig om van elk paar mee te nemen en niet enkel sample1 read1&2
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 



// ================================= ASSEMBLY  =============================
/* 
- Date creation: 03/04/2020
- Last adjustment: 3/6/2020
- Goal: Make de novo assembly with megahit and QA of the assembly with quast
- Faced problems:
    *problems with matplotlib (quast)
    *output assembly = final.contigs.fa, but then there's no samplename anymore!

Possibilities: include other assemblers and give a choice to users
*/

if(params.SE){
    //perform assembly with megahit for SE data
    process megahitSE {
        publishDir "$megahitdir", mode:'copy', overwrite: true
        
        input:
        file(trimfq) from trimmedSE_reads
        val(seq) from seq3a

        output:
        file("assembly/*.fa") into megahitSE 
        //path("assembly/*.fa") into megahitSE-chew 
        val(seq) into seq4a
        
        when:
        params.SE

        script:
        """
        megahit -r "${trimfq}" -o assembly --out-prefix ${trimfq.simpleName}-${seq} -t $params.cpu
        """
    }
}
else{
    //perform assembly with megahit for PE data
    process megahitPE {
        publishDir "$megahitdir", mode:'copy', overwrite: true
        
        input:
        set sample_id, file(trimfq) from trimmedPE_reads
        val(seq) from seq3b

        output:
        file("assembly/*.fa") into megahitPE  
        //path("assembly/*.fa") into megahitPE-chew  
        val(seq) into seq4b
        
        when:
        params.PE

        script:
        """
        megahit -1 "${trimfq[0]}" -2 "${trimfq[1]}" -o assembly --out-prefix ${sample_id}-${seq} -t $params.cpu
        """
    }
}
// split the file and stop division in SE and PE for quality control and further analysis with chewbbaca
if (params.SE){
    megahitSE.into {megahit_quast; megahit_mlst ;megahit_chew}
    //megahit_quast = megahitSE
    //megahit_chew = megahitSE-chew
    seq4a.into {seq5;seq5a}

}
else {
    megahitPE.into {megahit_quast; megahit_mlst; megahit_chew}
    //megahit_quast = megahitPE
    //megahit_chew = megahitPE-chew
    seq4b.into {seq5;seq5b}
}
//import matpotlib for quast
process matplotlib {
    script:
    """
    #!/usr/bin/env python3
    import matplotlib
    """
}
// Control of assembly with Quast: TO TEST
process quastSE{
    publishDir "$quastdir", mode:'copy', overwrite: true
       
    input:
    file(assembly) from megahit_quast
    val(seq) from seq5

    output:
    file("quast-${assembly.simpleName}/report.html")
    file("quast-${assembly.simpleName}/report.pdf")
    val(seq) into seq6

    script:
    """
    metaquast.py --threads "${task.cpus}" -o quast-${assembly.simpleName}-${seq} -s "${assembly}"
    """
}



// ================================= Scaffolding =============================
/* 
- Date creation: 03/04/2020
- Last adjustment: 3/6/2020
- Goal: 

TODOOOOOOOOOOOOOOOOOO
*/





// ================================= MLST type =============================
/* 
- Date creation: 26/06/2020
- Last adjustment: 26/06/2020
- Goal: Define ST-types for simple MLST analysis based on available schemes
*/

//Check MLST-type   TO TEST
process mlsttype{
        publishDir "$mlstdir", mode:'copy', overwrite: true

        input:
        file(assembly) from megahit_mlst

        output:
        file("*.tsv") 

        script:
        """
        mlst $megahitdir/assembly/*.fa > ST-types.tsv
        """
}



'''
// ================================= cgMLST analysis =============================
/* 
- Date creation: 3/6/2020
- Last adjustment: 3/6/2020
- Goal: Do actual cg/wgMLST analysis based on given scheme
- Faced problems
    * run stops because previous files are detected, however because nextflow is not interactive we 
      can not type to continue by "yes" -> fix by defining -fr = forced reset
    * analysis file by file but should provide a path to include ALL files at once
        ONly need to do this once!! -> Fix make input4 folder
    * make sure processes are done in sequential manner with maxForks 1
*/


//define the folder with the assemblies
//copy generated assemblies to the input folder for analysis
process input4a{

        input:
        file(assembly) from megahit_chew 
        val(seq) from seq6

        output:
        val(seq) into seq7
 
        script:
        """
        mkdir -p $outputdir
        mkdir -p $outputdir/Analysis
        mkdir -p $outputdir/Analysis/input4analysis
        cp ${assembly} $input4dir
        """
    }
//also copy given assemblies to the input folder for analysis
if (params.assem){
    Channel
    .fromPath("${assemdir}")
    .set {assem_ch}
    process input4b{

        input:
        file(assem) from assem_ch
        val(seq) from seq7

        output:
        val(seq) into seq8

        script:
        """
        cp ${assem} $input4dir
        """
    }
}


///TO FIX: problem with missing ANALYSIS FOLDER!!!


//do analysis
process analysis{

    publishDir "$analysisdir", mode:'copy', overwrite: true

    input:
    val(seq) from seq8

    output:
    file("*")
    file('*.tsv') into toclean

    script:
    """
    chewBBACA.py AlleleCall -i $input4dir -g $inputscheme -o Results --cpu $params.cpu --ptf $trainingfile --fr
    """
}

//cleanup  -- TO TEST
process cleanup{
    publishDir "$analysisdir", mode:'copy', overwrite: true

    input:
    result from toclean

    output:
    file("*")

    script:
    """
    chewBBACA.py ExtractCgMLST -i "${result}" -o Ready4MST
    """
}


//show some feedback about the pipeline  -- TO TEST
process feedback{
    echo true

    input:
    file(reportWG) from wgMLSTlogging
    file(counts) from ch_counts
    val(txid) from mtxid8

    script:
    """     
    #!/usr/bin/env python3
        
    print(" ")
    print ("REFERENCES USED FOR ANALYSIS")
    allcall= "$reportWG"
    txt1 = open (allcall, 'r')
    counter1 = 0
    for line1 in txt1:
        if counter1 == 3:
            print(line1, end='')
        counter1+=1
    txt1.close()
    print("out of ${params.count} requested")
    print(" ")
    """
}
'''

// =============================  Finished message ====================================
workflow.onComplete {
	log.info ( workflow.success ? "\nFINISHED!\n" : "Oops something went wrong!" )
}



//TO IMPROVE WHEN THERE'S EXCESS TIME
// Include non-zipped fastq's (check the given files and adjust protocol accordingly)
// Adjust nextflow.config
//  use docker container -> to set + make nextflow.config file!!

'''
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''