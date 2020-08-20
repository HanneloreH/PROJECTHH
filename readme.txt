
=====================================================================================
======================================= OUTB8 =======================================
=====================================================================================

Establishing a diagnostic pipeline for cgMLST analysis of outbreak samples

Author		:Hannelore Hamerlinck
Contact		:hannelore.hamerlinck@hotmail.com
Date		:13/08/2020



------------------------------------ INTRODUCTION ------------------------------------
Core genome multilocus sequence typing (cgMLST) is an elaborated method for characterizing isolates of microbial species on a DNA level using Whole Genome Sequencing (WGS). The high resolution typing of isolates has proven to be very useful for outbreak analysis. This project aims to automate the labour intensive analysis by generating pipelines that perform:
-	Generation of a cgMLST scheme based on NCBI reference sequences available for the chose organism
-	Performing cgMLST analysis on raw sequencing data (fastq) of suspected outbreak samples.

This OUTB8 project was executed as part of an internship "banaba bio-informatics" from October 2019 until August 2020.



---------------------------------- FOLDER STRUCTURE ---------------------------------
GITHUB FOLDER:
https://github.com/HanneloreH/PROJECTHH

CODE:
--	= Folder
---	= File
(! * !) = Most data in this folder is NOT shown on github (only locally) because of github filesize and memory restrictions

GENERAL STRUCTURE:
-- Data 		: (! * !) this folder contains all sorts of freely available and generated data 
   --  cgMLSTschemes 	: (! * !) this folder contains generated cgMLST schemes per txid
   --  dataKP-assembly	: (! * !) this folder contains assemblies from databases, projects...  
			from Klebsiella pneumoniae species
   --  RawData-KP	: (! * !) this folder contains raw fastq files from databases, projects...  
			from Klebsiella pneumoniae species
   --  refDB		: (! * !) this folder contains reference assemblies per given txid
   --  TrainingFiles	: this folder contains training files generated with prodigal

-- Notebook		: this folder contains all notes kept in the jupyter notebook
   --  Archive		: this folder contains archived notebook sheets
   --  DATA-Getting as..: this file states how the assemblies were obtained
   --  DATA-Getting Ra..: this file states how the fastq files were obtained
   --  FINAL TESTS	: this file gives a summary of the final pipeline tests on the local computer
   --  GIT-PROJECTHH..	: this file states how the github was made
   --  HPC		: this file gives the preparations of running hte pipeline on the HPC supercomputer
   --  PIPE Nextfl..	: this file states how the pipeline software (Nextflow) was installed
   --  PLOT-minimum..	: this file states how the script for the minimum spanning tree was created in R
   --  TOOLS-Chew..test.: this file states how the chewBBACA pipeline was tested with new data
   --  TOOLS-Chew..tut..: this file states how chewBBACA was installed and tested following a tutorial
   --  TOOLS install..  : this file states how software was installed and test (except chewBBACA)
REMARK: for notes on the pipeline itself please check Scripts/Pipeline/OUTB8-scheme.nf and Scripts/Pipeline/OUTB8-analysis.nf

-- Scripts		: this folder contains all generated scripts/pipelines...
   --  archive		: this folder contains archived scripts and information
   --  Input		: this folder contains small input files for testing (single/paired fastq & assemblies)
   --  Pipeline		: => => => this folder contains the generated pipelines and R script  <= <= <= 
   --  tests4scripts	: this folder contains scripts that were used for testing

-- Tests		: (! * !) this folder contains input and output of test-runs
   --  BAIT8-analysis-PE: "  *remark: BAIT8 = typo -> should have been OUTB8
   --  BAIT8-analysis-SE: "  *remark: BAIT8 = typo -> should have been OUTB8
   --  OUTB8-analysis-SE: "
   --  OUTB8-analysis-O.: "
   --  OUTB8-analysis-X.: "
   --  OUTB8-scheme-ex..: "
   --  OUTB8-scheme-fa..: "
   --  OUTB8-scheme-tx..: "

-- Tools		: (! * !) this folder contains instructions, information, code or test files for specific tools (chewBBACA, kraken, medusa, megahit, prodigal, quast, spades)

--- .gitignore		: this file contains all folders and files that were not uploaded to github
---  readme.txt		: readme-file




--------------------------------------- HOW TO --------------------------------------
1. ANALYSIS PIPELINE "OUTB8-analysis" v1.0
Do WGS bacterial analysis on fastq files (format *.fastq.gz) based on a known scheme for cg/wgMLST (includes trimming and assembly)

	- install following software:
	  * Nextflow v20.01.0.5264
	  * chewBBACA v2.5.4
	  * python3 v3.7.7
	  * fastqc v0.11.9
	  * multiqc v1.8
	  * fastp v0.20.0
	  * megahit v1.2.9
	  * matplotlib v3.3.0
	  * metaquast v5.0.2
	  * mlst v2.19.0
	  * R v3.6.3
	  * matplotlib v3.3.0
	- Download  OUTB8-analysis.nf and MST.R and run the script
		INPUT (must define):
		* --PE or --SE  :paired or single end data
		* --reads       :path to reads (if not in same folder)
		* --scheme      :path to scheme (folder with fasta files defined by scheme)
		OUTPUT: 
		* QA results before and after trimming in folder "Quality" 
		* simple MLST type analysis in folder "MLSTtypes"
		* Assembly and assembly quality parameters in folder "Assemblies-megahit"
		* Results and Minimum Spanning Tree in folder "Analysis"
		EXAMPLE:
		   $ nextflow run OUTB8-analysis.nf --PE --reads fastqs/ --scheme CGscheme/
		// Running requires prodigal training file (.trn) in current directory
	- optional parameters:
	    	--assem     OPTIONAL give path to established assemblies that need to be included in the analysis (format: *.fa)
	   	--cpu       give maximal number of CPUs (default = 1)
		--help      show help message
		--meta      give filename (.csv) of metadata, first column are the sample names IN CORRECT ORDER (= one line per sample), all following columns will be used as color-code in MST-plots
		--output    give path to output folder
		--training  give path to training file 
		--v         show version
		--x         if "true" pipeline will run starting from given assemblies (=cgMLST analysis only)
		-resume     use to continue an analysis that was run (partly) before



2. SCHEME PIPELINE "OUTB8-scheme" v1.0
Create a cg/wgMLST scheme for your WGS bacterial analysis based on NCBI reference assemblies based on txid OR 1 or more fastq samples (zipped)
	- install following software:
	  * Nextflow v20.01.0.5264
	  * chewBBACA v2.5.4
	  * python3 v3.7.7
	  * kraken2 v2.0.8-beta
	  * esearch v13.3
 	  * bit-dl-ncbi-assemblies
	  * prodigal V2.6.3

	- Download OUTB8-scheme.nf and run the script
		INPUT (must define):
		* --fastq    : give 1 or more fastq files (zipped)
		    OR
		* --txid     : give a known txid number (NCBI): https://www.ncbi.nlm.nih.gov/taxonomy
		OUTPUT: 
		* wgMLST scheme 
		* cgMLST scheme (default 95%)
		* prodigal training file
		EXAMPLE:
		   $ nextflow run OUTB8-scheme.nf --txid 573   
		   $ nextflow run OUTB8-scheme.nf --fastq "DNA123.fastq.gz"
	- optional parameters:
		 --count     give a maximum number of reference assemblies to use in the scheme (default 100)
				indicate "all" if you want to use all complete reference assemblies available at NCBI
		 --cpu       give maximal number of cpu's
		 --help      show help message
		 --krakendb  give path to kraken database 
		 --krakenout give path to folder where output of kraken2 can be stored (default current directory)
		 --perc      give percentage for cgMLST creation (default 0.95)
		 --refDB     give path to place/search reference assemblies (=output folder for reference assemblies)
		 --scheme    give path to wg/cg schemes (=output folder for the scheme)
		 --training  give path to training file (=output folder for the training file)
 		 --v         show version
		 -resume     use to continue an analysis that was run (partly) before

