#!/usr/bin/env python3

'''
def concat (file1, file2):
    import csv
    with open(file1, mode='r') as fil1:
        with open(file2, mode='r') as fil2:
           # with open ("concatenate.tsv", mode:'w') as concat:
                reader1= csv.reader(fil1, delimiter='\t')
                reader2= csv.reader(fil2, delimiter='\t')
                line_count = 1
                for row1 in reader1:
                    for col1 in row1:  #aantal kolommen definieren
                        for row2 in reader2:
                            for col2 in row2:

                        if (line_count<5):
                            print(col1)
                        line_count += 1
           # concat.close()
        fil2.close()
    fil1.close()
concat("/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/cgMLST_completegenomes/cgMLST.tsv", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/results/results_20200427T171943/results_alleles.tsv")
'''


#alternative: use pandas


def concat (file1, file2):
    import pandas
    fil1 = pandas.read_csv(file1, sep='\t', header=0)
    fil2 = pandas.read_csv(file2, sep='\t', header=0)
    concat = pandas.merge(fil2, on ='column_name')
    #write file

concat("/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/cgMLST_completegenomes/cgMLST.tsv", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/results/results_20200427T171943/results_alleles.tsv")