
=====================================================================================
======================================= OUTB8 =======================================
=====================================================================================

Establishing a diagnostic pipeline for cgMLST analysis of outbreak samples

Author		:Hannelore Hamerlinck
Contact		:hannelore.hamerlinck@hotmail.com
Date		:22/07/2020



------------------------------------ INTRODUCTION ------------------------------------
Core genome multilocus sequence typing (cgMLST) is an elaborated method for characterizing isolates of microbial species on a DNA level using Whole Genome Sequencing (WGS). The high resolution typing of isolates has proven to be very useful for outbreak analysis. This project aims to automate the labour intensive analysis by generating pipelines that perform:
-	Generation of a cgMLST scheme based on NCBI reference sequences available for the chose organism
-	Performing cgMLST analysis on raw sequencing data (fastq) of suspected outbreak samples.
- 	**under construction** Comparing different assemblers

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
   --- DATA-Getting as..: this file states how the assemblies were obtained
   --- DATA-Getting Ra..: this file states how the fastq files were obtained
   --- GIT-PROJECTHH..	: this file states how the github was made
   --- PIPE Nextfl..	: this file states how the pipeline software (Nextflow) was installed
   --- TOOLS-Chew..test.: this file states how the chewBBACA pipeline was tested with new data
   --- TOOLS-Chew..tut..: this file states how chewBBACA was installed and tested following a tutorial
   --- TOOLS install..  : this file states how software was installed and test (except chewBBACA)
   --- General info..	: this file contains information from before starting the pipeline generation (outdated)

-- Scripts		: this folder contains all generated scripts/pipelines...
   --  archive		: this folder contains archived scripts and information
   --  Input		: this folder contains small input files for testing (single/paired fastq & assemblies)
   --  Pipeline		: this folder contains the generated pipelines and R script
   --  tests4scripts	: this folder contains scripts that were used for testing

-- Tools		: (! * !) this folder contains instructions, information, code or test files for specific tools (chewBBACA, kraken, medusa, megahit, prodigal, quast, spades)

--- .gitignore		: this file contains all folders and files that were not uploaded to github
---  readme.txt		: readme-file



--------------------------------------- HOW TO --------------------------------------
1. ANALYSIS PIPELINE
Do WGS bacterial analysis on fastq files (format *.fastq.gz) based on a known scheme for cg/wgMLST (includes trimming and assembly)

	- install following software:
	  * Nextflow
	  * chewBBACA (install in desired python environment)
	  * fastqc
	  * multiqc
	  * fastp
	  * megahit
	  * matplotlib
	  * metaquast
	  * mlst
	  * Rscript
	- Download the OUTB8-analysis.nf and MSTfile and run the script
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
	   	--env       give path to python3 environment (running e.g. chewBBACA, matplotlib...)
		--help      show help message
		--meta      give filename of metadata
		--output    give path to output folder
		--training  give path to training file 
		--x         if "true" pipeline will run starting from given assemblies (=cgMLST analysis only)
		-resume     use to continue an analysis that was run (partly) before






