#!/bin/bash
source activate /data/gent/vo/001/gvo00121/env/envOUTB8
nextflow run OUTB8-analysis-4HPC.nf --PE --reads /data/gent/vo/001/gvo00121/nextflow-test/Test2-input/input/ --scheme /data/gent/vo/001/gvo00121/nextflow-test/Test2-input/scheme --output /data/gent/vo/001/gvo00121/nextflow-test/Test2-input/output --cpu 8 --meta /data/gent/vo/001/gvo00121/nextflow-test/Test2-input/input/meta.csv

