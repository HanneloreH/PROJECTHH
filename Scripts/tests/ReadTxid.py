#!/usr/bin/env python3

#read txid from text file
import csv


with open('report-test.txt', mode='r') as report:
    txtreader = csv.reader(report, delimiter='\t')
    line_count = 1
    for row in txtreader:
        if(row[3]=="S"):
            if(line_count<=1):           
                print(row[4])
                line_count += 1
report.close()
