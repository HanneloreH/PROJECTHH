#!/usr/bin/python3

#isntall eutils: pip3 install --user eutils AND pip install eutils
from eutils import Client

eclient = Client(api_key="f7bba41c57271397bbc91985bce382bf6***")
print("\nUsing NCBI E-utilities in Python\n")

'''#does work
gene_esearch = eclient.esearch(db='gene',term='TNF')
#does not work
assembly_esearch = eclient.esearch(db='assembly',term='klebsiella pneumoniae')

'''

esearch -db assembly -query '573[txid] AND "complete genome"[filter] AND "latest refseq"[filter]'