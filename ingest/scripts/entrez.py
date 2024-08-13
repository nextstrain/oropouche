from Bio import Entrez
from sys import stdin,stdout,stderr
import pandas as pd

Entrez.email = "hello@nextstrain.org"

if __name__=="__main__":
    accessions = list(pd.read_csv(stdin,sep='\t',index_col='accession').index.values)
    print(f"{len(accessions) = }. Querying entrez...", file=stderr)
    # accessions=accessions[:10] # useful for testing
    handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text")
    for line in handle.readlines():
        print(line, end='', file=stdout)
    handle.close()
