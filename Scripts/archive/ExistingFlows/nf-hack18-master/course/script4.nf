/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

println """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/* 
 * create a transcriptome file object given then transcriptome string parameter
 */
transcriptome_file = file(params.transcriptome)
 
/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    file transcriptome from transcriptome_file
     
    output:
    file 'index' into index_ch

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

// PE reads matchen en kijken of het empty is dan opslaan als channel
Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Oops! Cannot find any file matching: ${params.reads}"  }
    .set { read_pairs_ch } 

//read_pairs_ch.println()


process quantification {
    publishDir "$baseDir/chunks", mode: 'copy', overwrite: false

    input:
    file index from index_ch
    set pair_id, file(reads) from read_pairs_ch
 
//second input is declared as a set composed by two elements: the pair_id and the reads in order to match the structure of the items emitted by the read_pairs_ch channel.


    output:
    file(pair_id) into quant_ch
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """

}

quant_ch.println()
