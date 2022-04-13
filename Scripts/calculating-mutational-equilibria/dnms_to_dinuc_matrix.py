### IMPORT STATEMENTS ###

import os
import random
import csv
from generic import run_in_parallel
import numpy as np

### SOURCE FILES ###

source = './DNMs/decode_DNMs.tsv'
chromsomes_source = './Chromosomes/WGS/Unzipped/Homo_sapiens.GRCh38.dna.chromosome.X.fa'
out_source = './CSV/dinucleotide_matrix.tsv'

### MAIN ###

def main():

    #Load DNMs
    raw_file = open(os.path.abspath(source)).read()
    rawSplit = raw_file.split('\n')[1:]

    #Filter to remove indels
    mutations = []
    for i in rawSplit[1:]:
        if i != '':
            line = i.split('\t')
            ref = line[2]
            alt = line[3]
            if len(ref) == 1 and len(alt) == 1:
                mutations.append(i.replace('chr', '').split('\t'))

    #Iterate through mutations and infer dinucleotide changes
    dicts = []

    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(mutations, ['foo'], function, workers=workers)

    total_success = 0
    total_fail = 0

    for process in processes:
        output, success, fail = process.get()
        dicts.append(output)
        total_success += success
        total_fail += fail

    total_changes_dict = {}
    for d in dicts:
        for k, v in d.items():
            if k in total_changes_dict:
                total_changes_dict[k] += v
            else:
                total_changes_dict[k] = v

    #Convert dictionary to nested list
    nested_results = [[k.split('_')[0].lower(), k.split('_')[1].lower(), str(v)] for k,v in total_changes_dict.items()]

    #Output to TSV file
    filepath = os.path.abspath(out_source)
    headers = ['ref', 'alt', 'count']

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(i for i in headers)
        for j in nested_results:
            writer.writerow(j)

### FUNCTIONS ###

def function(mutations):

    changes = {}
    success = 0
    fail = 0

    for mutation in mutations:

        chr = mutation[0]
        location = int(mutation[1])
        ref = mutation[2]
        alt = mutation[3]

        #Load the relevant chromosome
        chromosome_file = open(os.path.abspath(chromsomes_source.replace('X', chr))).read()
        sequence = "".join(chromosome_file.split('\n')[1:])
        nucleotide_check = sequence[location-1]
        if ref == nucleotide_check:
            success+= 1
            dinucleotide_1_ref = sequence[location-2] + ref
            dinucleotide_1_alt = sequence[location-2] + alt
            change1 = dinucleotide_1_ref + '_' + dinucleotide_1_alt
            dinucleotide_2_ref = ref + sequence[location]
            dinucleotide_2_alt = alt + sequence[location]
            change2 = dinucleotide_2_ref + '_' + dinucleotide_2_alt

            #Record dinucleotide change 1
            if change1 not in changes:
                changes[change1] = 1
            else:
                changes[change1] += 1

            #Record dinucleotide change 2
            if change2 not in changes:
                changes[change2] = 1
            else:
                changes[change2] += 1

            print (dinucleotide_1_ref, 'to', dinucleotide_1_alt, '+', dinucleotide_2_ref, 'to', dinucleotide_2_alt)
        else:
            fail += 1

    return changes, success, fail

### RUN ###

if __name__ == '__main__':
    main()
