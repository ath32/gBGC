### IMPORT STATEMENTS ###

import os
import re
import numpy as np
import csv
import random
from tqdm import tqdm
from generic import run_in_parallel

### SOURCE FILES ###

full_source = './Primates/Ortho/Aligned'
high_source = './Primates/Ortho/Aligned_GC_rich'
low_source = './Primates/Ortho/Aligned_GC_poor'

### MAIN ###

def main():

    #Get genomic flux rates and pTGA
    genomic_taa_tga, genomic_tga_taa, genomic_pTGA = get_nulls(full_source)

    #Simulate GC rich set - see methods for methodology
    high_sims = []
    workers = int(os.cpu_count()) - 1
    sims = list(range(1, 101))
    high_processes = run_in_parallel(sims, ['foo', high_source, genomic_taa_tga, genomic_tga_taa], get_sims, workers=workers)
    for process in high_processes:
        output = process.get()
        high_sims.extend(output)

    #Simulate GC poor set - see methods for methodology
    low_sims = []
    workers = int(os.cpu_count()) - 1
    sims = list(range(1, 101))
    low_processes = run_in_parallel(sims, ['foo', low_source, genomic_taa_tga, genomic_tga_taa], get_sims, workers=workers)
    for process in low_processes:
        output = process.get()
        low_sims.extend(output)

    #Select random pairs and calculate the difference in pTGA, compare to the real sets
    counter = 0
    no_sims = 0
    for i in tqdm(range(1,10001)):
        random_high = np.random.choice(high_sims)
        random_low = np.random.choice(low_sims)
        diff = random_high - random_low
        observed = 0.08544058 #Edit this observed difference according to the trio tested, sorry this is hard coded
        if diff > observed:
            counter += 1
        no_sims += 1

    #P = the number with a higher pTGA difference than the real gene sets divided by number of sims
    p = counter / no_sims
    print (p)

### FUNCTIONS ###

def get_nulls(full_source):

    #Generate genomic pTGA
    taa = 0
    tga = 0
    taa_tga = 0
    tga_taa = 0

    for root, dirs, filenames in os.walk(full_source):
        for f in tqdm(filenames):
            if f.endswith(".fa"):

                #Get aligned file (fasta)
                path = os.path.join(full_source, f)
                raw_file = open(path).read()

                #Get list of sequences (third one is the outgroup)
                seq_list = []
                split = raw_file.split('>')
                for i in split[1:]:
                    sections = i.split('\n')
                    seq_sections = sections[1:]
                    seq = ''.join(seq_sections)
                    seq_list.append(seq.upper())

                #Count switches by position
                in1 = seq_list[0]
                in2 = seq_list[1]
                out = seq_list[2]

                #Get genes that we can infer an ancestral state
                in1_codons = [in1[i:i+3] for i in range(0, len(in1), 3)]
                in2_codons = [in2[i:i+3] for i in range(0, len(in2), 3)]
                out_codons = [out[i:i+3] for i in range(0, len(out), 3)]

                stops = ['TAA', 'TGA', 'TAG']

                #Find position of the ancestral stop codon from the alignment and get the end branch states
                for i,n in enumerate(out_codons):
                    if n in stops and in1_codons[i] in stops and in2_codons[i] in stops:
                        if out_codons[i-1] and out_codons[i-2] not in stops:
                            if in1_codons[i-1] and in1_codons[i-2] not in stops:
                                if in2_codons[i-1] and in2_codons[i-2] not in stops:

                                    codon_index = i

                                    branch1 = in1_codons[i]
                                    branch2 = in2_codons[i]

                                    #Now get the ancestral reconstructed stop
                                    try:
                                        a_path = os.path.join(full_source, f + '.state')
                                        a_file = open(a_path).read()
                                        split = a_file.split('p_T\n')[1].split('\n')

                                        start_index = codon_index * 3
                                        end_index = start_index + 3

                                        rel_lines = split[start_index:end_index]
                                        reconstructed = "".join([i.split('\t')[2] for i in rel_lines])

                                        if reconstructed == 'TAA':
                                            taa += 1
                                            if branch1 == 'TGA' or branch2 == 'TGA':
                                                taa_tga += 1

                                        elif reconstructed == 'TGA':
                                            tga += 1
                                            if branch1 == 'TAA' or branch2 == 'TAA':
                                                tga_taa += 1

                                    except:
                                        continue

    genomic_taa_tga = taa_tga / taa
    genomic_tga_taa = tga_taa / tga
    genomic_pTGA = 1 / (1 + (genomic_tga_taa / genomic_taa_tga))

    return genomic_taa_tga, genomic_tga_taa, genomic_pTGA

def get_sims(sims, source, genomic_taa_tga, genomic_tga_taa):

    #According to genomic flux rates, simulate pTGA scores
    results = []
    for sim in sims:

        taa = 0
        tga = 0
        taa_tga = 0
        tga_taa = 0

        for root, dirs, filenames in os.walk(source):
            for f in tqdm(filenames):
                if f.endswith(".fa"):

                    #Get aligned file (fasta)
                    path = os.path.join(source, f)
                    raw_file = open(path).read()
                    gene_id = raw_file.replace('.fa', '')

                    #Get list of sequences (third one is the outgroup)
                    seq_list = []
                    split = raw_file.split('>')
                    for i in split[1:]:
                        sections = i.split('\n')
                        seq_sections = sections[1:]
                        seq = ''.join(seq_sections)
                        seq_list.append(seq.upper())

                    ids = []
                    for i in split:
                        if i != '':
                            splt = i.split('\n')
                            id = splt[0].split(';')[0].replace('ID=gene:', '')
                            ids.append(id)

                    #Count switches by position
                    in1 = seq_list[0]
                    in2 = seq_list[1]
                    out = seq_list[2]

                    #Get genes that we can infer an ancestral state
                    in1_codons = [in1[i:i+3] for i in range(0, len(in1), 3)]
                    in2_codons = [in2[i:i+3] for i in range(0, len(in2), 3)]
                    out_codons = [out[i:i+3] for i in range(0, len(out), 3)]

                    stops = ['TAA', 'TGA', 'TAG']

                    #Find position of the ancestral stop codon from the alignment and get the end branch states
                    for i,n in enumerate(out_codons):
                        if n in stops and in1_codons[i] in stops and in2_codons[i] in stops:
                            if out_codons[i-1] and out_codons[i-2] not in stops:
                                if in1_codons[i-1] and in1_codons[i-2] not in stops:
                                    if in2_codons[i-1] and in2_codons[i-2] not in stops:

                                        codon_index = i

                                        branch1 = in1_codons[i]
                                        branch2 = in2_codons[i]

                                        #Now get the ancestral reconstructed stop
                                        for root, dirs, filenames in os.walk(source):
                                            for f1 in filenames:
                                                if f1.endswith(".state"):
                                                    if f.split('.')[0] in f1:

                                                        #Get aligned file (fasta)
                                                        path = os.path.join(source, f1)
                                                        raw_file = open(path).read()
                                                        split = raw_file.split('p_T\n')[1].split('\n')

                                                        start_index = codon_index * 3
                                                        end_index = start_index + 3

                                                        rel_lines = split[start_index:end_index]
                                                        reconstructed = "".join([i.split('\t')[2] for i in rel_lines])

                                                        if reconstructed == 'TAA':
                                                            taa += 1
                                                            random_number = np.random.choice(['A', 'B'], replace=True, p=[genomic_taa_tga, 1-genomic_taa_tga])
                                                            if random_number == 'A':
                                                                taa_tga += 1

                                                        if reconstructed == 'TGA':
                                                            tga += 1
                                                            random_number = np.random.choice(['A', 'B'], replace=True, p=[genomic_tga_taa, 1-genomic_tga_taa])
                                                            if random_number == 'A':
                                                                tga_taa += 1

        null_taa_tga = taa_tga / taa
        null_tga_taa = tga_taa / tga
        null_pTGA =  1 / (1 + (null_tga_taa / null_taa_tga))
        results.append(null_pTGA)

    return results

### RUN ###

if __name__ == '__main__':
    main()
