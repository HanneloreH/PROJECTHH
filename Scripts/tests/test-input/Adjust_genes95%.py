#!/usr/bin/env python3


def Addpath (file, path):
    import csv

    with open(file, mode='r') as genes:
        txtreader = csv.reader(genes, delimiter='\t')
        line_count = 1
        for row in txtreader:
            print(row[4])
            line_count += 1
    genes.close()

Addpath(Genes_95%.txt, "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/schema/")