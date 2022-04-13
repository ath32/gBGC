### IMPORT STATEMENTS ###

import os
import random
import csv
from generic import run_in_parallel
import numpy as np

### SOURCE FILES ###

source = './DNMs/decode_DNMs_filtered_singles_gc.tsv'
chromsomes_source = './Chromosomes/WGS/Unzipped/Homo_sapiens.GRCh38.dna.chromosome.X.fa'
out_source = './CSV/mononucleotide_matrix.tsv'

### MAIN ###

def main():

    #load DNMs
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

    #Now iterate through mutations and infer mononucleotide changes, also record total nucleotide counts
    changes = []
    dinucleotides = []

    total_a = 0
    total_c = 0
    total_g = 0
    total_t = 0

    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(mutations, ['foo'], function, workers=workers)

    for process in processes:
        output, total_dinucs, a, c, g, t = process.get()
        changes.append(output)
        dinucleotides.append(total_dinucs)
        total_a += a
        total_c += c
        total_g += g
        total_t += t

    #Create dictionary containing the number of each type of mutation
    total_changes_dict = {}
    for d in changes:
        for k, v in d.items():
            if k in total_changes_dict:
                total_changes_dict[k] += v
            else:
                total_changes_dict[k] = v

    #Convert dictionary to nested list
    nested_results = [[k.split('_')[0].lower(), k.split('_')[1].lower(), str(v)] for k,v in total_changes_dict.items()]

    #Output file
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
    total_dinucs = {}
    a = 0
    c = 0
    g = 0
    t = 0

    for mutation in mutations:

        chr = mutation[0]
        location = int(mutation[1])
        ref = mutation[2]
        alt = mutation[3]
        gc = float(mutation[4])*100

        #load the relevant chromosome
        chromosome_file = open(os.path.abspath(chromsomes_source.replace('X', chr))).read()
        sequence = "".join(chromosome_file.split('\n')[1:])
        nucleotide_check = sequence[location-1]
        if ref == nucleotide_check and ref != 'N' and alt != 'N':

            change = ref + '_' + alt

            #Get total counts
            sequence = sequence[location-5001: location+4999]
            dinucleotides = []
            for i in range(len(sequence)-1):
                dinucleotide = sequence[i:i+2]
                dinucleotides.append(dinucleotide)
            a += sequence.count('A')
            c += sequence.count('C')
            g += sequence.count('G')
            t += sequence.count('T')

            for dinuc in dinucleotides:
                if dinuc not in total_dinucs:
                    total_dinucs[dinuc] = 1
                else:
                    total_dinucs[dinuc] += 1

            #Record change 1
            if change not in changes:
                changes[change] = 1
            else:
                changes[change] += 1

            print (ref, 'to', alt)

        else:
            print ('REF does not match hg19')

    return changes, total_dinucs, a, c, g, t

### RUN ###

if __name__ == '__main__':
    main()
