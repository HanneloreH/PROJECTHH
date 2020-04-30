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
#download via pip wa snot suffcient, also done via conda install pandas!

def concat (file1, file2):
    import pandas as pd #does work in terminal
    fil1 = pd.read_csv(file1, sep='\t', header=0)
    fil2 = pd.read_csv(file2, sep='\t', header=0)
    files = [fil1, fil2]
    combine = pd.concat(files, join="inner", ignore_index=True)
    #combineD = combine.drop(combine.columns[0], axis=1)
    #if deleting column 0 = FILE, how to delete -1?
    combine.to_csv("output.tsv", sep='\t', index=False) #added index=False

    

concat("/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/cgMLST_completegenomes/cgMLST.tsv", "/home/hannelore/PROJECTHH/Tools/chewBBACA_run1/results/results_20200427T171943/results_alleles.tsv")