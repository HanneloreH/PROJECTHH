 #!/usr/bin/env python3
 
def MLSTcopy (filo,inputpath,outputpath):
    import os
    import shutil
    fileList = os.listdir(inputpath)
    cfil = open (filo, 'r')
    countCG=0
    for line in cfil:
        countCG +=1
        line=line[:-1]
        countWG = 0
        for item in fileList:
            countWG +=1
            if(line == item):
                print(line)
        counts = open("counts.txt", "a")
        print()
    counts.write("number of loci: {}".format(countWG))
    counts.write("\nnumber of loci: {}".format(countCG))
    counts.close()
    cfil.close()

tekst= '/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-287-all/cgMLST/cgMLSTschema.txt'
inputtie = '/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-287-all/wgMLST/schema-287-all/'
outputtie ='/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-287-all/cgMLST/scheme-287-all-cgMLST/'

MLSTcopy(tekst, inputtie, outputtie)