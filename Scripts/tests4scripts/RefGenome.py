#!/usr/bin/python3



# Reading data from UCSC
print("\nREADING DATA FROM UCSC DATABASE")
db = my.connect(host="genome-euro-mysql.soe.ucsc.edu",
   user="genomep",
   passwd="password",
   db="hg19")
c = db.cursor()
no_rows = c.execute("""SELECT * FROM ensGene""")
# Fetch one row (ENST + ENSG)
print(c.fetchone())
print("Total rows ensGene: {}".format(no_rows))
# Total rows ensGene: 204940
