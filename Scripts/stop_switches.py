''' Stop switch analysis for the primary stop codon site. '''


### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict
import csv

### CHOOSE SOURCE FOLDER - Unhash the older of interest ###

#ALL
source = 'Arabidopsis/Ortho/Aligned_10'

def main():

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

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

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

            #Filtering criteria this time... no gaps in 2 upstream codons
            for i,n in enumerate(out_codons):
                if n in stops and in1_codons[i] in stops and in2_codons[i] in stops:
                    if out_codons[i-1] and out_codons[i-2] not in stops:
                        if in1_codons[i-1] and in1_codons[i-2] not in stops:
                            if in2_codons[i-1] and in2_codons[i-2] not in stops:
                                all_count += 1
                                #TGA genes
                                if n == 'TGA':
                                    tga_n += 1
                                    if in1_codons[i] == 'TAA' or in2_codons[i] == 'TAA':
                                        tga_taa += 1
                                    if in1_codons[i] == 'TAG' or in2_codons[i] == 'TAG':
                                        tga_tag += 1
                                #TAA genes
                                if n == 'TAA':
                                    taa_n += 1
                                    if in1_codons[i] == 'TGA' or in2_codons[i] == 'TGA':
                                        taa_tga += 1
                                    if in1_codons[i] == 'TAG' or in2_codons[i] == 'TAG':
                                        taa_tag += 1
                                #TAG genes
                                if n == 'TAG':
                                    tag_n += 1
                                    if in1_codons[i] == 'TAA' or in2_codons[i] == 'TAA':
                                        tag_taa += 1
                                    if in1_codons[i] == 'TGA' or in2_codons[i] == 'TGA':
                                        tag_tga += 1

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

    print (taa_tga_r + taa_tag_r + tga_taa_r + tga_tag_r + tag_taa_r + tag_tga_r)

    # Get null frequencies

    #Count switches
    aa_ga = 0
    aa_ag = 0
    ga_aa = 0
    ga_ag = 0
    ag_aa = 0
    ag_ga = 0

    ta_tg = 0
    tg_ta = 0

    #n counts
    aa_n = 0
    ag_n = 0
    ga_n = 0
    ta_n = 0
    tg_n = 0

    #First, get synonymous sites
    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get aligned file (fasta)
            path = os.path.join(source, f)
            raw_file = open(path).read()

            #Get list of sequences (third one is the outgroup)
            sequence_list = []
            split = raw_file.split('>')
            for i in split[1:]:
                sections = i.split('\n')
                seq_sections = sections[1:]
                seq = ''.join(seq_sections)
                seq_coding = list(seq[:18].upper())
                seq_syn = seq_coding[::3]
                sequence_list.append(seq.upper())

            #Count switches by position
            in1 = sequence_list[0]
            in2 = sequence_list[1]
            out = sequence_list[2]

            #Get genes that we can infer an ancestral state
            in1_codons = [in1[i:i+3] for i in range(0, len(in1), 3)]
            in2_codons = [in2[i:i+3] for i in range(0, len(in2), 3)]
            out_codons = [out[i:i+3] for i in range(0, len(out), 3)]

            stops = ['TAA', 'TGA', 'TAG']

            #Filtering criteria this time... no gaps in 2 upstream codons
            for i,n in enumerate(out_codons):
                if n in stops and in1_codons[i] in stops and in2_codons[i] in stops:

                    #Get 3' codons (after stop)
                    downstream_out = out_codons[i+1:]
                    downstream_in1 = in1_codons[i+1:]
                    downstream_in2 = in2_codons[i+1:]
                    #Turn into nucleotides
                    nts_out = "".join(downstream_out)
                    nts_in1 = "".join(downstream_in1)
                    nts_in2 = "".join(downstream_in2)

                    #Get all possible dinucleotide combinations
                    in1_dinuc = []
                    for i in range(len(nts_in1)-1):
                        dinuc = nts_in1[i:i+2]
                        in1_dinuc.append(dinuc)

                    in2_dinuc = []
                    for i in range(len(nts_in2)-1):
                        dinuc = nts_in2[i:i+2]
                        in2_dinuc.append(dinuc)

                    out_dinuc = []
                    for i in range(len(nts_out)-1):
                        dinuc = nts_out[i:i+2]
                        out_dinuc.append(dinuc)

                    #Get switch counts
                    for i,n in enumerate(out_dinuc):

                        #AA -> GA or AG
                        if n == 'AA':
                            if in1_dinuc[i] == 'AA' or in2_dinuc[i] == 'AA':
                                aa_n += 1
                                if in1_dinuc[i] == 'AA' and in2_dinuc[i] == 'AG':
                                    aa_ag += 1
                                if in1_dinuc[i] == 'AG' and in2_dinuc[i] == 'AA':
                                    aa_ag += 1
                                if in1_dinuc[i] == 'AA' and in2_dinuc[i] == 'GA':
                                    aa_ga += 1
                                if in1_dinuc[i] == 'GA' and in2_dinuc[i] == 'AA':
                                    aa_ga += 1
                        #AG -> AA or GA
                        if n == 'AG':
                            if in1_dinuc[i] == 'AG' or in2_dinuc[i] == 'AG':
                                ag_n += 1
                                if in1_dinuc[i] == 'AG' and in2_dinuc[i] == 'AA':
                                    ag_aa += 1
                                if in1_dinuc[i] == 'AA' and in2_dinuc[i] == 'AG':
                                    ag_aa += 1
                                if in1_dinuc[i] == 'AG' and in2_dinuc[i] == 'GA':
                                    ag_ga += 1
                                if in1_dinuc[i] == 'GA' and in2_dinuc[i] == 'AG':
                                    ag_ga += 1
                        #GA -> AA or AG
                        if n == 'GA':
                            if in1_dinuc[i] == 'GA' or in2_dinuc[i] == 'GA':
                                ga_n += 1
                                if in1_dinuc[i] == 'GA' and in2_dinuc[i] == 'AA':
                                    ga_aa += 1
                                if in1_dinuc[i] == 'AA' and in2_dinuc[i] == 'GA':
                                    ga_aa += 1
                                if in1_dinuc[i] == 'GA' and in2_dinuc[i] == 'AG':
                                    ga_ag += 1
                                if in1_dinuc[i] == 'AG' and in2_dinuc[i] == 'GA':
                                    ga_ag += 1

                        #TA -> TG
                        if n == 'TA':
                            if in1_dinuc[i] == 'TA' or in2_dinuc[i] == 'TA':
                                ta_n += 1
                                if in1_dinuc[i] == 'TA' and in2_dinuc[i] == 'TG':
                                    ta_tg += 1
                                if in1_dinuc[i] == 'TG' and in2_dinuc[i] == 'TA':
                                    ta_tg += 1
                        #TG -> TA
                        if n == 'TG':
                            if in1_dinuc[i] == 'TG' or in2_dinuc[i] == 'TG':
                                tg_n += 1
                                if in1_dinuc[i] == 'TG' and in2_dinuc[i] == 'TA':
                                    tg_ta += 1
                                if in1_dinuc[i] == 'TA' and in2_dinuc[i] == 'TG':
                                    tg_ta += 1
    #Frequencies
    aa_ga_f = aa_ga / aa_n
    aa_ag_f = aa_ag / aa_n
    ga_aa_f = ga_aa / ga_n
    ga_ag_f = ga_ag / ga_n
    ag_aa_f = ag_aa / ag_n
    ag_ga_f = ag_ga / ag_n
    ta_tg_f = ta_tg / ta_n
    tg_ta_f = tg_ta / tg_n

    aa_ga_r = aa_ga_f / (aa_ga_f + aa_ag_f + ga_aa_f + ga_ag_f + ag_aa_f + ag_ga_f)
    aa_ag_r = aa_ag_f / (aa_ga_f + aa_ag_f + ga_aa_f + ga_ag_f + ag_aa_f + ag_ga_f)
    ga_aa_r = ga_aa_f / (aa_ga_f + aa_ag_f + ga_aa_f + ga_ag_f + ag_aa_f + ag_ga_f)
    ga_ag_r = ga_ag_f / (aa_ga_f + aa_ag_f + ga_aa_f + ga_ag_f + ag_aa_f + ag_ga_f)
    ag_aa_r = ag_aa_f / (aa_ga_f + aa_ag_f + ga_aa_f + ga_ag_f + ag_aa_f + ag_ga_f)
    ag_ga_r = ag_ga_f / (aa_ga_f + aa_ag_f + ga_aa_f + ga_ag_f + ag_aa_f + ag_ga_f)

    #Write CSV
    csv_total = []
    csv_total.append(['Switch', 'Type', 'Freq', 'Hits', 'Total'])
    csv_total.append(['TAA > TGA', 'Stop switch', taa_tga_r, taa_tga, taa_n])
    csv_total.append(['TAA > TGA', 'Null', aa_ga_r, aa_ga, aa_n])
    csv_total.append(['TAA > TAG', 'Stop switch', taa_tag_r, taa_tag, taa_n])
    csv_total.append(['TAA > TAG', 'Null', aa_ag_r, aa_ag, aa_n])
    csv_total.append(['TGA > TAA', 'Stop switch', tga_taa_r, tga_taa, tga_n])
    csv_total.append(['TGA > TAA', 'Null', ga_aa_r, ga_aa, ga_n])
    csv_total.append(['TGA > TAG', 'Stop switch', tga_tag_r, tga_tag, tga_n])
    csv_total.append(['TGA > TAG', 'Null', ga_ag_r, ga_ag, ga_n])
    csv_total.append(['TAG > TAA', 'Stop switch', tag_taa_r, tag_taa, tag_n])
    csv_total.append(['TAG > TAA', 'Null', ag_aa_r, ag_aa, ag_n])
    csv_total.append(['TAG > TGA', 'Stop switch', tag_tga_r, tag_tga, tag_n])
    csv_total.append(['TAG > TGA', 'Null', ag_ga_r, ag_ga, ag_n])

    filename = "Switches_10_relative.csv"
    subdir = "Arabidopsis/CSV"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
