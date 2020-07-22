

params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir= "$baseDir/output"

println """reads: $params.reads 
transcriptome: $params.transcriptome
multiqc: $params.multiqc
Output directory: $params.outdir
"""
