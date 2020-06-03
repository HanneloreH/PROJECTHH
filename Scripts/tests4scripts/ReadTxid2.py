#!/usr/bin/env python3
#read txid from text file

def GetTxid(report):
    import csv
    with open(report, mode='r') as rap:
        txtreader = csv.reader(rap, delimiter='\t')
        line_count = 1
        for row in txtreader:
            if(row[3]=="S"):
                if(line_count<=1):           
                    print(row[4])
                    line_count += 1
    rap.close()
