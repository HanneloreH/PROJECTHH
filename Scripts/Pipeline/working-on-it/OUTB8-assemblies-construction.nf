#!/usr/bin/env nextflow
 
/*
========================================================================================
                         cgMLST analysis on fastq data
========================================================================================
SUMMARY
Do several assemblies on trimmed fastq data

INPUT (must define):
    * --PE or --SE  :paired or single end data
    * --reads       :path to trimmed fastq files

OUTPUT: 
    * Assembly for every sample and assembly quality parameters in folder "Assemblies-*"

EXAMPLE INPUT:
    * nextflow run Pipeline/PROJECTHH-assemblies.nf --PE --reads "/home/hannelore/PROJECTHH/Data/RawData-KP/" 


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
params.cpu= 1
params.help = false
params.output= "$baseDir/output"
params.PE = false
params.reads = "$baseDir/*{1,2}.fastq.gz"
params.SE = false


// Define all folders 
//Remove final "/" from workingfolders and define them:
//function to trim final /
def trimFolder = { 
    it.endsWith("/") ? it[0..-2] : it
}
//input folders
fastq_f = trimFolder("$params.reads")
fastq_files = fastq_f + "/*{1,2}.fastq.gz"

//output folders
outputdir = trimFolder("$params.output")
megahitdir = outputdir + "/Assemblies-megahit"
spadesdir = outputdir + "/Assemblies-spades"
velvetdir = outputdir + "/Assemblies-velvet"

//resultdir = outputdir + "/Results"  //defined in analysis

// Transform reads to files and form pairs and set channels

//Input channel for further analysis (split between PE and SE reads with if-clause)
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

//split files in 3 for 3 assemblies
 read_pairs_ch.into {read_pairsMEGA_ch; read_pairsSPADES_ch; read_pairsVELVET_ch}


// Print the parameters set
println """\
        Analysis starting using following parameters:   
        ==============================================================
        * Number of CPUs                  : ${params.cpu}
        * Output-folder                   : ${outputdir}
        * Folder with input reads         : ${fastq_f}
        * PE reads                        : ${params.PE}
        * SE reads                        : ${params.SE}
         """
         .stripIndent()



// ================================= ASSEMBLY  =============================
/* 
Make assembly from megahit, spades and velvet
*/


//Megahit
if(params.SE){
    //perform assembly with megahit for SE data
    process megahitSE {
        publishDir "$megahitdir", mode:'copy', overwrite: true
        
        input:
        file(trimfq) from read_pairsMEGA_ch

        output:
        file("assembly/*.fa")
        
        when:
        params.SE

        script:
        """
        megahit -r "${trimfq}" -o assembly --out-prefix ${trimfq.simpleName} -t $params.cpu
        """
    }
}
else{
    //perform assembly with megahit for PE data
    process megahitPE {
        publishDir "$megahitdir", mode:'copy', overwrite: true
        
        input:
        set sample_id, file(trimfq) from read_pairsMEGA_ch

        output:
        file("assembly/*.fa")
        
        when:
        params.PE

        script:
        """
        megahit -1 "${trimfq[0]}" -2 "${trimfq[1]}" -o assembly --out-prefix ${sample_id} -t $params.cpu
        """
    }
}


//Spades
if(params.SE){
    //perform assembly with spades for SE data
    process spadesSE {
        publishDir "$spadesdir", mode:'copy', overwrite: true
        
        input:
        file(trimfq) from read_pairsSPADES_ch

        output:
        file("assembly-spades-${trimfq}/*")
        
        when:
        params.SE

        script:
        """
        spades.py --se ${trimfq}  -o assembly-spades-${trimfq} -t $params.cpu
        """
    }
}
else{
    //perform assembly with spades for PE data
    process spadesPE {
        publishDir "$spadesdir", mode:'copy', overwrite: true
        
        input:
        set sample_id, file(trimfq) from read_pairsSPADES_ch

        output:
        file("assembly-spades-${sample_id}/*")
        
        when:
        params.PE

        script:
        """
        spades.py -1 "${trimfq[0]}" -2 "${trimfq[1]}" -o assembly-spades-${sample_id} -t $params.cpu
        """
    }
}

//Velvet
if(params.SE){
    //perform assembly with velvet for SE data
    process velvetSE {
        publishDir "$velvetdir", mode:'copy', overwrite: true
        
        input:
        file(trimfq) from read_pairsVELVET_ch

        output:
        file("assembly-velvet-${trimfq}/*")
        
        when:
        params.SE

        script:
        """
        velveth velvet21-${trimfq} 21 -fastq.gz -short ${trimfq}
        velvetg
        """
    }
}
else{
    //perform assembly with velvet for PE data
    process velvetPE {
        publishDir "$velvetdir", mode:'copy', overwrite: true
        
        input:
        set sample_id, file(trimfq) from read_pairsMEGA_ch

        output:
        file("assembly-velvet-${sample_id}/*")
        
        when:
        params.PE

        script:
        """
        velveth velvet21-${sample_id} 21 -fastq.gz -shortPaired ${trimfq[0]} ${trimfq[1]}
        """
    }
}






// =============================  Finished message ====================================
workflow.onComplete {
	log.info ( workflow.success ? "\nFINISHED!\n" : "Oops something went wrong!" )
}


// TODO: 
// find better way to reach adapter file
// test for zipped or not fastqs
// nextflow.config
//  use docker container (?)-> to set + make nextflow.config file!!
// database with already downloaded assemblies per txid
'''
do -with-docker
ALTERNATIVE in nexflow.config file:
docker.enabled = true
'''