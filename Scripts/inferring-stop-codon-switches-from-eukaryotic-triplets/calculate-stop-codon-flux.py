### IMPORT STATEMENTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict
import csv

### SOURCE FILES ###

source = 'Primates/Ortho/Aligned_GC_rich'

### MAIN ###

def main():

    #Counters for switches
    taa_tga = 0
    taa_tag = 0
    tga_taa = 0
    tga_tag = 0
    tag_taa = 0
    tag_tga = 0

    #n counts
    taa_n = 0
    tga_n = 0
    tag_n = 0

    #Infer stop codon switches
    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            if f.endswith(".fa"):

                #Get aligned file (fasta)
                path = os.path.join(source, f)
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
                                                        taa_n += 1
                                                        if branch1 == 'TGA' or branch2 == 'TGA':
                                                            taa_tga += 1
                                                        elif branch1 == 'TAG' or branch2 == 'TAG':
                                                            taa_tag += 1

                                                    elif reconstructed == 'TGA':
                                                        tga_n += 1
                                                        if branch1 == 'TAA' or branch2 == 'TAA':
                                                            tga_taa += 1
                                                        elif branch1 == 'TAG' or branch2 == 'TAG':
                                                            tga_tag += 1

                                                    elif reconstructed == 'TAG':
                                                        tag_n += 1
                                                        if branch1 == 'TAA' or branch2 == 'TAA':
                                                            tag_taa += 1
                                                        elif branch1 == 'TGA' or branch2 == 'TGA':
                                                            tag_tga += 1

    #From counts, calculate frequencies
    taa_tga_f = (taa_tga / taa_n)
    taa_tag_f = (taa_tag / taa_n)
    tga_taa_f = (tga_taa / tga_n)
    tga_tag_f = (tga_tag / tga_n)
    tag_taa_f = (tag_taa / tag_n)
    tag_tga_f = (tag_tga / tag_n)

    #From frequencies, calculate relative frequencies
    taa_tga_r = taa_tga_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    taa_tag_r = taa_tag_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tga_taa_r = tga_taa_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tga_tag_r = tga_tag_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tag_taa_r = tag_taa_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tag_tga_r = tag_tga_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)

    #Write CSV, nulls are redundant
    csv_total = []
    csv_total.append(['Switch', 'Type', 'Freq', 'Hits', 'Total'])
    csv_total.append(['TAA > TGA', 'Stop switch', taa_tga_f, taa_tga, taa_n])
    csv_total.append(['TAA > TGA', 'Null', 'N/a', 'N/a', 'N/a'])
    csv_total.append(['TAA > TAG', 'Stop switch', taa_tag_f, taa_tag, taa_n])
    csv_total.append(['TAA > TAG', 'Null', 'N/a', 'N/a', 'N/a'])
    csv_total.append(['TGA > TAA', 'Stop switch', tga_taa_f, tga_taa, tga_n])
    csv_total.append(['TGA > TAA', 'Null', 'N/a', 'N/a', 'N/a'])
    csv_total.append(['TGA > TAG', 'Stop switch', tga_tag_f, tga_tag, tga_n])
    csv_total.append(['TGA > TAG', 'Null', 'N/a', 'N/a', 'N/a'])
    csv_total.append(['TAG > TAA', 'Stop switch', tag_taa_f, tag_taa, tag_n])
    csv_total.append(['TAG > TAA', 'Null', 'N/a', 'N/a', 'N/a'])
    csv_total.append(['TAG > TGA', 'Stop switch', tag_tga_f, tag_tga, tag_n])
    csv_total.append(['TAG > TGA', 'Null', 'N/a', 'N/a', 'N/a'])

    #Write output file - name according to GC-rich/GC-poor or HRG/LRG
    filename = "GC-rich-switches.csv"
    subdir = "./Primates/CSV"
    filepath = os.path.join(subdir, filename)
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
