
#!/usr/bin/env python3

import os
import shutil

inputtie = '/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-548/wgMLST/schema-548/'
outputtie = '/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-548/cgMLST/scheme-548-cgMLST/'
fileList = os.listdir(inputtie)
print(fileList)

cfil = open ('/home/hannelore/PROJECTHH/Data/cgMLSTschemes/MLST-548/cgMLST/cgMLSTschema.txt', 'r')
for line in cfil:
    #Had problems matching line to file: but found that the number of characters differed (1 extra for line)
    #tried to remove with strip and replace (" ","") but this did not work, removing the last character did, it was an /n
    """
    line.split()
    line.replace(" ","")
    """
    line=line[:-1]
    for item in fileList:
        '''
        print("line is", line)
        print(type(line))
        print(len(line))
        print("item is", item)
        print(type(item))
        print(len(item))
        '''
        if(line == item):
            inputpath= str(inputtie) + str(line)
            print(inputpath)
            outputpath= str(outputtie) + str(line)
            print(outputpath)
            shutil.copyfile(str(inputpath), str(outputpath))
cfil.close()
