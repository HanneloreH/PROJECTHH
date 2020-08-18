#!/bin/bash
source activate /data/gent/vo/001/gvo00121/envOUTB8
nextflow run OUTB8-analysis-tryPREP.nf --SE --reads ~/PROJECTHH/Tests/BAIT8-analysis-SE/input2/ --output ~/PROJECTHH/Tests/BAIT8-analysis-SE/output-prep --scheme ~/PROJECTHH/Data/cgMLSTschemes/MLST-573-c20-Prepped/cgMLST/scheme-573-c20-cgMLST-prep --cpu 2  --training ~/PROJECTHH/Data/TrainingFiles/txid573-c20.trn --meta ~/PROJECTHH/Tests/BAIT8-analysis-SE/input/meta.csv

