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
                #lineB= line2 + "_short.fasta_bsr.txt"

                inputtoA= str(inputpath) + "short/" + str(lineA)
                outputtoA= str(outputpath) + "short/" + str(lineA)
                shutil.copyfile(str(inputtoA), str(outputtoA))

                #inputtoB= str(inputpath) + "short/" + str(lineB)
                #outputtoB= str(outputpath) + "short/" + str(lineB)
                #shutil.copyfile(str(inputtoB), str(outputtoB))
    counts = open("counts.txt", "a")
    counts.write("number of loci: {}".format(countCG))
    counts.close()            
    cfil.close()
    

tekst= '/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-573-c25/wgMLST/Genes_95%.txt'
inputtie = '/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-573-c25/wgMLST/schema-573-c25/'
outputtie ='/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-573-c25/cgMLST/scheme-573-c25-cgMLST/'


MLSTcopy(tekst, inputtie, outputtie)
