### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict
import csv
import scipy.stats as st
from tqdm import tqdm

### CHOOSE SOURCE FOLDER - Unhash the older of interest ###

#ALL
source = 'Primates/Ortho/Aligned_R2'

def main():

    files = []
    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            if f.endswith(".fa"):
                print ('adding file...')
                files.append(f)

    print ('calculating observed...')
    ptga_o = get_switches(files)
    print (ptga_o, 'Ob. done')

    #Bootstrap!
    sims = 0

    ptga_list = []

    print ('doing sims...')

    while sims < 1000:
        sample = np.random.choice(files,len(files),replace=True)
        tga_e = get_switches(sample)
        ptga_list.append(tga_e)
        sims += 1
        print (sims, 'done')

    ptga_ci = st.t.interval(alpha=0.95, df=len(ptga_list)-1, loc=np.mean(ptga_list), scale=st.sem(ptga_list))
    ptga_se = np.std(ptga_list, ddof=1) / np.sqrt(np.size(ptga_list))
    print (ptga_o, ptga_se, ptga_ci)

def get_switches(list):

    #Gene counts
    all_count = 0

    #Count switches
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

    for f in list:
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
                                try:
                                    a_path = os.path.join(source, f + '.state')
                                    a_file = open(a_path).read()
                                    split = a_file.split('p_T\n')[1].split('\n')

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

                                except:
                                    continue

    taa_tga_f = (taa_tga / taa_n)
    taa_tag_f = (taa_tag / taa_n)
    tga_taa_f = (tga_taa / tga_n)
    tga_tag_f = (tga_tag / tga_n)
    tag_taa_f = (tag_taa / tag_n)
    tag_tga_f = (tag_tga / tag_n)

    taa_tga_r = taa_tga_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    taa_tag_r = taa_tag_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tga_taa_r = tga_taa_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tga_tag_r = tga_tag_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tag_taa_r = tag_taa_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)
    tag_tga_r = tag_tga_f / (taa_tga_f + taa_tag_f + tga_taa_f + tga_tag_f + tag_taa_f + tag_tga_f)

    ptga =  1 / (1 + (tga_taa_r / taa_tga_r))

    return ptga


### RUN ###

if __name__ == '__main__':
    main()
