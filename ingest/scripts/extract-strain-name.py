from sys import stdin,stdout,stderr
from Bio import SeqIO
from collections import defaultdict, Counter

def extract(feature, key):
    if key in feature.qualifiers:
         return feature.qualifiers[key][0]
    return False

if __name__=="__main__":
    n_records = 0
    n_strain_names_found = 0
    names = Counter()

    print("accession\tstrain", file=stdout)
    for record in SeqIO.parse(stdin, "genbank"):
        n_records+=1
        accession = record.id.split('.')[0]
        source = [feature for feature in record.features if feature.type == 'source'][0]
        strain_name = extract(source, 'strain') or extract(source, 'isolate')
        if strain_name:
            n_strain_names_found+=1
            strain_name = strain_name.replace(" ","_") # replace spaces in strain names with underscores
            print(f"{accession}\t{strain_name}", file=stdout)
            names.update([strain_name])
        else:
            print(f"{accession}\t{accession}", file=stdout)

    print(f"Found strain names in {n_strain_names_found}/{n_records} records after inspecting 'strain' or 'isolate' fields. Using 'accession' for the rest.", file=stderr)

    print("Collecting matching strain names together resulted in the following counts:", file=stderr)
    occurances = defaultdict(int)
    for c in names.values():
        occurances[c]+=1
    for n in sorted(occurances.keys()):
        print(f"\t{n=} {occurances[n]} times" ,file=stderr)
