#!/usr/bin/env python3

'''
#TRY1 FAILED: could not find right delimiter (not " ", not "\t", not ",")

def Addpath (file, path):
    import csv

    with open(file, mode='r') as genes:
        txtreader = csv.reader(genes, delimiter=' ')
        line_count = 1
        for row in txtreader:
            if(line_count<=5):
                print(row)
                line_count += 1
    genes.close()

Addpath("Genes_95%.txt", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/schema/")
'''
'''
#TRY2
def Addpath (file, path):
    fil = open (file, 'r')
    output = open("Genes_95%-PATH.txt", 'w+')
    line_count = 1
    for line in fil:
        for word in line.split():
            if(line_count>3):
                output.write("{}{}\n".format(path,word))
            line_count += 1

    output.close()
    fil.close()

Addpath("Genes_95%.txt", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/schema/")
'''


'''
#TRY3: only word
def Addpath (file, path):
    fil = open (file, 'r')
    output = open("Genes_95%-SEP.txt", 'w+')
    line_count = 1
    for line in fil:
        for word in line.split():
            if(line_count>3):
                output.write("{}\n".format(word))
            line_count += 1

    output.close()
    fil.close()

Addpath("Genes_95%.txt", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/schema/")'''


#TRY4: word + schema
def Addpath (file, path):
    fil = open (file, 'r')
    output = open("Genes_95%-path.txt", 'w+')
    line_count = 1
    for line in fil:
        for word in line.split():
            if(line_count>3):
                output.write("schema/{}\n".format(word))
            line_count += 1

    output.close()
    fil.close()

Addpath("Genes_95%.txt", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/schema/")