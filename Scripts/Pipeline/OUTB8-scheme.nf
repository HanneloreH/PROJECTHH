#!/usr/bin/env nextflow
 
/*
========================================================================================
                OUTB8-scheme v1.0 : cgMLST scheme creation pipeline
========================================================================================
SUMMARY
Create a cg/wgMLST scheme for your WGS bacterial analysis by using a number of NCBI reference assemblies based on a txid OR a fastq file

INPUT (must define):
    *--fastq    : give 1 or more fastq files (zipped)
    OR
    *--txid     : give a known txid number (NCBI): https://www.ncbi.nlm.nih.gov/taxonomy

OUTPUT: 
    *wgMLST scheme 
    *cgMLST scheme (default 95%)
    *prodigal training file

EXAMPLE INPUT:
    *nextflow run OUTB8-scheme.nf --txid 573   
    *nextflow run OUTB8-scheme.nf --fastq "DNA123.fastq.gz"

AUTHOR
Hannelore Hamerlinck <hannelore.hamerlinck@hotmail.com>

----------------------------------------------------------------------------------------
*/



// ====================================  Show help message and version =========================================
/* 
- Goal: Show the help message when asked for help, give all possible parametesr
- Date creation: 15/05/2020
- Last adjusted: 23/07/2020
*/

def helpMessage() {
    log.info"""
    ==================== OUTB8-scheme ====================

    Create a cg/wgMLST scheme for your WGS bacterial analysis 
    based on NCBI reference assemblies based on txid OR 1 or more fastq samples (zipped)

    example1:   nextflow run OUTB8-scheme.nf --txid "817"
    example2:   nextflow run OUTB8-scheme.nf --fastq "DNA123.fastq.gz"

    --count     give a maximum number of reference assemblies to use in the scheme (default 100)
                indicate "all" if you want to use all complete reference assemblies available at NCBI
    --cpu       give maximal number of cpu's
    --env       give path to python3 environment (running e.g. chewBBACA...)
    --fastq     give one or more fastq files (zipped)
    --help      show help message
    --krakendb  give path to kraken database 
    --krakenout give path to folder where output of kraken2 can be stored (default current directory)
    --perc      give percentage for cgMLST creation (default 0.95)
    --refDB     give path to place/search reference assemblies (=output folder for reference assemblies)
    --scheme    give path to wg/cg schemes (=output folder for the scheme)
    --training  give path to training file (=output folder for the training file)
    --txid      give a known txid number (NCBI): https://www.ncbi.nlm.nih.gov/taxonomy
    --v         show version
    -resume     use to continue an analysis that was run (partly) before

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

def version(){
    println("OUTB8-scheme v1.0")
}

if (params.v){
    version()
    exit 0
}


// ============================================= Set Parameters ========================================
/* 
- Goal: define all default parameters and folders (can by adjusted by user in command prompt)
- Date creation: 15/05/2020
- Last adjusted: 26/06/2020
- Faced problems:
    * problem cutting of / at the end (first script took the first / -> unknown folder)
      Fix: use [0..-2]
*/
 
// Define general parameters (defaults)
params.count = 100
params.cpu = 1  //$task.cpus does not work...
params.fastq = "zero"
params.help = false
params.krakendb = "/home/hannelore/PROJECTHH/Tools/kraken-db/16S_SILVA138_k2db" //SILVA, TO ADJUST!
params.krakenout = "$baseDir/outputKraken"
params.perc = 0.95
params.refDB = "/home/hannelore/PROJECTHH/Data/refDB"
params.scheme = "/home/hannelore/PROJECTHH/Data/cgMLSTschemes"
params.training= "/home/hannelore/PROJECTHH/Data/TrainingFiles"
params.v = false
params.txid = "zero"

//Remove final "/" from workingfolders and define them:
//function to trim final /
def trimFolder = { 
    it.endsWith("/") ? it[0..-2] : it
}
//do it for working folders
fastqx = trimFolder("$params.fastq")
trainingx = trimFolder("$params.training")
refDBx = trimFolder("$params.refDB")
cgMLSTx = trimFolder("$params.scheme")
refDBa = refDBx + "/accession-lists"
krakenx = trimFolder("$params.krakenout")
krakenxres = krakenx + "/reports"
krakenxmini = krakenx + "/mini-fastqs"
krakendbx = trimFolder("$params.krakendb")
envdir = trimFolder("$params.env")
envactivate = envdir + "/bin/activate"

 // ===========================  Set Python3 environment ===========================
 
process environment {
    script:
    """
    source $envactivate
    """
}




// ===============================================  Get txid ===========================================
/* 
- Goal: Divide the flow based on the user choice for --txid or --fastq to get the correct txid 
  When a fastq is given kraken2 analysis is performed on a fraction of the reads (10000 = mini-fastq) to determine the txid.
- Remark: Created TXID for use as txid value, create txida for further analysis
- Date creation: 15/05/2020 (* processes were already written in april'20)
- Last adjusted: 26/06/2020
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
    *Extracting txid to use as value from channel (after kraken2)
        Tried many things, some examples below;
            * tried to save the txid as file and get text: TXID2 = file("printtxid.txt").getText('UTF-8')
            * save txid as view: TXID =  gtxid1.view {it.trim()}
            * use subscribe: TXID = gtxid1.subscribe { println "value: $it" } -> println "please please please: ${TXID}"
            * using collethttps://www.nextflow.io/docs/latest/operator.html?highlight=first#transforming-operators: 
            ...
        -> I was not able to fix this, so I gave message to restart the pipeline with the correct txid, is needed anyhow if there are multiple TXIDs
*/


if (params.txid=="zero"){
    //=check if no txid is given

    if (params.fastq!="zero"){
        //= check if fastq is given

        //Set fastq file
        fastq_files = fastqx + "/*.fastq.gz"
        Channel
             .fromPath( fastq_files)
             .ifEmpty { error "INPUT ERROR: Cannot find any fastq.gz files, please try again"} //=double warning...
             .set { fastq_ch }
        
         //Make miniFastq (10000 reads)
        process reduce4kraken {

            input:
            file(fastqFile) from fastq_ch
            output:
            file("${fastqFile}.mini.fastq.gz") into reduced4kraken
            when:
            params.fastq
            script:
            """
            zcat ${fastqFile} | head -40000 | gzip > ${fastqFile}.mini.fastq.gz
            """
        }
        
        //Perform kraken2 to get txid from mini-fastq
        process kraken2 {
            publishDir "$krakenxmini", mode:'copy', overwrite: true, pattern: "${miniFQ}"
            publishDir "$krakenxres", mode:'copy', overwrite: true, pattern: "Kraken2-report*"
            
            input:
            file(miniFQ) from reduced4kraken

            output:
            file("Kraken2-report-${miniFQ.simpleName}.txt") into report
            file("${miniFQ}") into filenames
            
            when:
            params.fastq

            script:
            """
            kraken2 --use-names --report Kraken2-report-${miniFQ.simpleName}.txt --threads $params.cpu \
            -db $krakendbx --gzip-compressed ${miniFQ}
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
            column -s, -t < !{rep}| awk '$4 == "G"'| head -n 1 | cut -f5 | tr -dc '0-9'
            '''
            //Now searches for "G" = Genus (because of silva database) -> to adjust to "S" = species
        }

        process printtxid {
            echo true
            publishDir "$krakenx", mode:'copy', overwrite: true

            input:
            val(txid) from gtxid
            file(name) from filenames

            output:
            file("*-TXID.txt") into finish

            script:
            """
            echo 'Kraken2 calculation on 10000 reads for ${name.simpleName} gives txid ${txid} as dominant species' | tee -a ${name.simpleName}-TXID.txt
            """
        }

        process finishing {
            echo true

            script:
            """
            echo "\n-------------- Please restart the pipeline with the corret TXID --------------"
            exit 0
            """
        }

    }

    else {
        //=if no fastq is given
        println "INPUT ERROR: please provide either a txid OR fastq file with --txid OR --fastq"
    }
}

else {
    //=if txid is given use this one, print txid for user
    TXID = params.txid
    println """\nThe given txid is ${TXID}, number of reference assemblies is MAXIMUM ${params.count}""".stripIndent()

    //write variable into channel to use in processes
    txida = Channel.value(params.txid)       
}


//only execute the next steps if TXID was defined by user
if (params.txid!="zero"){


    // =================================== Get Reference Assemblies ====================================
    /*
    - Goal: Get reference assemblies from NCBI via txid, 
    all proceses in this section wil run if there is no txid sepcfic accession.txt file available in the DB-folder (if-clause)
    - Date creation: 18/05/2020 (* processes were already written in april'20)
    - Last adjusted: 25/05/2020
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
            fix: NEED $TXID (see previous)
    */

    //define suffix:
    if ("${params.count}"=="all"){
        suffix = "${params.count}"
    }
    else{
        suffix = "c${params.count}"
    }


    //check if the accession exists with given count
    accFile = file("${refDBa}/${TXID}accessions-c${params.count}.txt")
    if (accFile.isFile() ) { 

        println "\nSEARCHING ASSEMBLIES... The accession file is already available continue to scheme generation "
        txidF1=txida
    }
    else {
        println "\nSEARCHING ASSEMBLIES... NO accession file available for the given count, start acquiring reference assemblies"

        // get full accession-list ("all") file from NCBI with esearch (only complete assemblies)
        process accessions {
            publishDir "$refDBa", mode:'copy', overwrite: true
                
            input:
            val txid from txida

            output:
            file("${txid}accessions-all.txt") into accessionlistALL
            val(txid) into acctxid

            script:
            """
            mkdir -p $refDBa
            esearch -db assembly -query '${txid}[txid] AND "complete genome"[filter] AND "latest refseq"[filter]' \
            | esummary | xtract -pattern DocumentSummary -element AssemblyAccession > ${txid}accessions-all.txt
            """
        }

        //limit number of assemblies to count ONLY if NOT indicated "all" at params.count
        if (params.count != "all") { 
            
            process limitAssemblies{
                
                publishDir "$refDBa", mode:'copy', overwrite: true
                    
                input:
                val(txid) from acctxid
                file(access) from accessionlistALL //not necessary anymore

                output:
                file("${txid}accessions-${suffix}.txt") into accessionlistCOUNT
                val(txid) into acctxid2

                """     
                #!/usr/bin/env python3

                txtfile= "$refDBa/${txid}accessions-all.txt"
                output = open("${txid}accessions-${suffix}.txt", "a")
                counter = 0
                max = ${params.count}
                txt = open (txtfile, 'r')

                for line in txt:
                    if counter < max:
                        output.write (line) 
                    counter +=1
                txt.close()
                """
            }
        }
        else{
        accessionlistCOUNT = accessionlistALL
        acctxid2 = acctxid
        }


        // get the fasta files matching the accession numbers and save in refDB
        process refAssembly {
            publishDir "$refDBx", mode:'copy', overwrite: true
            
            input:
            file(accesions) from accessionlistCOUNT
            val(txid) from acctxid2

            output:
            file("refDB-${txid}-${suffix}/GC*")
            val(txid) into acctxid3

            //when:
            //when is not necessary because it has already been defined

            script:
            """
            mkdir -p refDB-${txid}-${suffix}
            cp ${accesions} refDB-${txid}-${suffix}
            cd refDB-${txid}-${suffix}
            bit-dl-ncbi-assemblies -w ${accesions} -f fasta -j 10
            """

        }

        // unzip the fasta files
        process unzip {
            publishDir "$refDBx/refDB-${txid}-${suffix}/", mode:'copy', overwrite: true
            
            input:
            val(txid) from acctxid3

            output:
            file("unzipped/GC*")
            val(txid) into acctxid4

            //when:
            //when is not necessary because it has already been defined

            script:
            """
            mkdir unzipped
            cp $refDBx/refDB-${txid}-${suffix}/*.gz unzipped/
            cd unzipped
            gunzip *.gz
            """
        }
        //XTRA: look into possiblity to NOT save unzipped files (to spare memory) only in temporary files!

        // make multifasta from ref assemblies
        process multi {
            publishDir "$trainingx", mode:'copy', overwrite: true
            
            input:
            val(txid) from acctxid4

            output:
            file("multi-txid${txid}-${suffix}.fa")
            val(txid) into txidF1

            //when:
            //when is not necessary because it has already been defined by accesions!

            script:
            """
            cat $refDBx/refDB-${txid}-${suffix}/unzipped/GCF* > multi-txid${txid}-${suffix}.fa 
            """
        }
    }



    // =================================== creating wgMLST ====================================
    /*
    - Goal: Create a wgMLST scheme based on loci from reference assemblies
    - Date creation: 19/05/2020
    - Last adjusted: 26/06/2020
    - Faced problems:
            * Problems trying to make processes run sequential but loose from each other (based on what's already available)
            FIX:treat cgMLST analysis as one entity, and stop the script if the scheme already exists
            * Problem: wgMSLT analysis and allele calling: time consuming step!
            * Problem: paralogs: cannot predict name of folder under allelecalling
            FIX: use "*"
            * Problem short.fasta_bsr.txt files are not all created during wgMLSt analysis, only part and of 1 genome
                - is necessary for allelecalling!
            * Problem: Allelecalling does not happen on all fragments only a small part, unclear why????
    */

    //end script if cg/wg-analysis is available in the given database
    cgFile = file("${cgMLSTx}/MLST-${TXID}-${suffix}/cgMLST/cgMLST.tsv")
    if (cgFile.isFile() ) {
        println """\nCREATING cgMLST SCHEME: The scheme for txid ${TXID}-${suffix} has already been created, stopping workflow""".stripIndent()
        exit 0
    }


    //create trainingfile with prodigal
    process training {
        publishDir "$trainingx", mode:'copy', overwrite: true
        
        input:
        val(txid) from txidF1

        output:
        file("txid${txid}-${suffix}.trn")
        val(txid) into mtxid1

        script:
        """
        prodigal -i $trainingx/multi-txid${txid}-${suffix}.fa -p single -t txid${txid}-${suffix}.trn 
        """
    }

    //create wgMLST scheme: 
    //REMARK: TIME CONSUMING (>> depending on number of genomes)
    process wgMLST {
        publishDir "$cgMLSTx/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid1

        output:
        file("MLST-${txid}-${suffix}-PRE/wgMLST/*")
        val(txid) into mtxid2a

        script:
        """
        mkdir MLST-${txid}-${suffix}-PRE
        cd MLST-${txid}-${suffix}-PRE
        mkdir wgMLST
        cd wgMLST
        chewBBACA.py CreateSchema -i $refDBx/refDB-${txid}-${suffix}/unzipped/ \
        -o schema-${txid}-${suffix}-PRE --cpu $params.cpu --ptf $trainingx/txid${txid}-${suffix}.trn
        """
    }

    //Added extra step: preparation of scheme, because of compatibility issues with chewBBACA v2.1.0 and v2.5.0
    process prepscheme {
        publishDir "$cgMLSTx/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid2a

        output:
        file("MLST-${txid}-${suffix}/wgMLST/*")
        val(txid) into mtxid2b

        script:
        """
        mkdir MLST-${txid}-${suffix}
        cd MLST-${txid}-${suffix}
        mkdir wgMLST
        cd wgMLST
        chewBBACA.py CreateSchema -i $cgMLSTx/MLST-${txid}-${suffix}-PRE/wgMLST/schema-${txid}-${suffix}-PRE/ \
        -o schema-${txid}-${suffix} --cpu $params.cpu --ptf $trainingx/txid${txid}-${suffix}.trn
        """
    }


    // allele calling 1 of reference assemblies on wgMLST scheme
    // REMARK: TIME CONSUMING (>> depending on number of genomes)
    // REMARK: IN PIPELINE GOES WRONG: only part of the loci are used:  the ones with short fasta bsr! -> problem lays in wgMLST
    process allelcall {
        publishDir "$cgMLSTx/MLST-${txid}-${suffix}/wgMLST/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid2b

        output:
        file("*")
        file("allelecalling/*/logging_info.txt") into wgMLSTlogging
        val(txid) into mtxid3

        script:
        """
        mkdir allelecalling
        chewBBACA.py AlleleCall -i $refDBx/refDB-${txid}-${suffix}/unzipped/ \
        -g $cgMLSTx/MLST-${txid}-${suffix}/wgMLST/schema-${txid}-${suffix}/ -o allelecalling/ \
        --cpu $params.cpu 
        """
        // training file was removed because of compatibility errors: "--ptf $trainingx/txid${txid}-${suffix}.trn"
    }

    // Remove paralogs
    process paralogs {
        publishDir "$cgMLSTx/MLST-${txid}-${suffix}/wgMLST/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid3

        output:
        file("alleleCallMatrix_cg.tsv")
        val(txid) into mtxid4

        script:
        """
        chewBBACA.py RemoveGenes -i $cgMLSTx/MLST-${txid}-${suffix}/wgMLST/allelecalling/*/results_alleles.tsv \
        -g $cgMLSTx/MLST-${txid}-${suffix}/wgMLST/schema-${txid}-${suffix}/allelecalling/*/RepeatedLoci.txt \
        -o alleleCallMatrix_cg
        """
    }

    // Check quality of wgMLST scheme (visual in html)
    process gquality {
        publishDir "$cgMLSTx/MLST-${txid}-${suffix}/wgMLST/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid4

        output:
        file("GenomeQualityPlot.html")
        file("Genes_95%.txt")
        file("removedGenomes.txt")
        val(txid) into mtxid5
        
        script:
        """
        chewBBACA.py TestGenomeQuality -i $cgMLSTx/MLST-${txid}-${suffix}/wgMLST/alleleCallMatrix_cg.tsv -n 13 -t 200 -s 5
        """
    }

    // Define loci present in params.perc % of reference genomes (default= 95%)
    //also prepare folders for next process (e.g. short)
    process cgMLST {
        publishDir "$cgMLSTx/MLST-${txid}-${suffix}/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid5

        output:
        file("cgMLST/*")
        val(txid) into mtxid6
        
        script:
        """
        mkdir cgMLST
        mkdir cgMLST/scheme-${txid}-${suffix}-cgMLST
        mkdir cgMLST/scheme-${txid}-${suffix}-cgMLST/short
        chewBBACA.py ExtractCgMLST -i $cgMLSTx/MLST-${txid}-${suffix}/wgMLST/alleleCallMatrix_cg.tsv \
        -o cgMLST -p $params.perc
        """
    }

    //Copy loci that are present in 95% of the reference genomes into cgMLST folder with python script 
    //+ also copy the corresponding "short-files"
    process cgMLSTloci{
        echo true

        input:
        val(txid) from mtxid6

        output:
        file("counts.txt") into ch_counts
        val(txid) into mtxid7
        
        """
        #!/usr/bin/env python3
        def MLSTcopy (filo,inputpath,outputpath):
            import os
            import shutil
            fileList = os.listdir(inputpath)
            cfil = open (filo, 'r')
            countCG=-1
            for line in cfil:
                countCG +=1
                line=line[:-1]
                for item in fileList:
                    if(line == item):
                        #move basis
                        inputto= str(inputpath) + str(line)
                        outputto= str(outputpath) + str(line)
                        shutil.copyfile(str(inputto), str(outputto))
                    
                        #move files in folder short (short was created in process cgMLST)
                        line2= line[:-6]
                        lineA= line2 + "_short.fasta"
                        lineB= line2 + "_short.fasta_bsr.txt"

                        inputtoA= str(inputpath) + "short/" + str(lineA)
                        outputtoA= str(outputpath) + "short/" + str(lineA)
                        shutil.copyfile(str(inputtoA), str(outputtoA))

                        inputtoB= str(inputpath) + "short/" + str(lineB)
                        outputtoB= str(outputpath) + "short/" + str(lineB)
                        shutil.copyfile(str(inputtoB), str(outputtoB))
            counts = open("counts.txt", "a")
            counts.write("number of loci: {}".format(countCG))
            counts.close()            
            cfil.close()

        tekst= '$cgMLSTx/MLST-${txid}-${suffix}/cgMLST/cgMLSTschema.txt'
        inputtie = '$cgMLSTx/MLST-${txid}-${suffix}/wgMLST/schema-${txid}-${suffix}/'
        outputtie ='$cgMLSTx/MLST-${txid}-${suffix}/cgMLST/scheme-${txid}-${suffix}-cgMLST/'

        MLSTcopy(tekst, inputtie, outputtie)
        """
    }

    //Evaluate the created scheme 
    process schemeEval {
        publishDir "$cgMLSTx/MLST-${txid}-${suffix}/", mode:'copy', overwrite: true
        
        input:
        val(txid) from mtxid7

        output:
        file("Scheme-Evaluation/*")
        val(txid) into mtxid8

        script:
        """
        mkdir "Scheme-Evaluation"
        chewBBACA.py SchemaEvaluator -i $cgMLSTx/MLST-${txid}-${suffix}/cgMLST/scheme-${txid}-${suffix}-cgMLST/ \
        -ta 11 -l Scheme-Evaluation/Evaluation-${txid}-${suffix}.html --cpu $params.cpu
        """
    }

    //show some feedback about the pipeline
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

        print ("WG-MLST SCHEME")
        allcall= "$reportWG"
        txt2 = open (allcall, 'r')
        counter2 = 0
        for line2 in txt2:
            if counter2 == 4:
                print(line2, end='')
            counter2+=1
        txt2.close()
        print(" ")

        print ("CG-MLST SCHEME")
        MLST= "$counts"
        txt3 = open (MLST, 'r')
        counter3 = 0
        for line3 in txt3:
            if counter3 == 0:
                print(line3, end='')
            counter3+=1
        txt3.close()
        print(" ")
        print(" ")

        """
    }


    // =============================  Finished message ====================================
    workflow.onComplete {
        log.info ( workflow.success ? "\ncgMLST SCHEME FOR ${TXID} CREATED SUCCESSFULLY!\n\n" : "Oops something went wrong, please check your input again!" )
    }


    //trial with txid 817: bact. fragilis (only 13 ref assemblies)


    //TO IMPROVE WHEN THERE'S EXCESS TIME
    // Software in container (DOCKER)
    // TXID variable from fastq file (as value) so scheme creation can start immediately
    // Adjust the flow to include existing pipelines (do conversion to chewBBACA + quality check)
    // Fix problem when interupting the pipeline and removing wgMLST (fastq-bsr files are not copied again)
    // Adjust "sufix" value to actual value (e.g. c13 in stead of c100 in txid 817)
    // Make file out of feedback
}


