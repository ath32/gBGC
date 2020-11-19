''' Take sequences, calculate observed and expected (dinuc. controlled simulations)
    for 'G' containing codons. '''


### IMPORTS ###

import os
import numpy as np
import csv
import time
import scipy.stats as st
from Bio.Seq import Seq

### CHOOSE SOURCE FOLDER (EMBL) ###

source = 'TAG_Tests/Test2/FASTA/downstream_data.txt'

### FUNCTIONS ###

def main():
    folder_parse(source)

def folder_parse(source):

    csv_total = []

    path = os.path.abspath(source)
    raw = open(path).read()

    #Start timer
    start_time = time.time()

    #Split into gene list for use in functions
    genes = raw.strip().split('>')
    genes_c = list(filter(None, genes))
    gene_list = []

    total_cds = ''
    for i in genes_c:
        sequence = "".join(i.split('\n')[1:]).upper()
        if sequence[0:3] == 'ATG':
            if sequence[len(sequence)-3:len(sequence)] == 'TAA' or sequence[len(sequence)-3:len(sequence)] == 'TGA' or sequence[len(sequence)-3:len(sequence)] == 'TAG':
                gene_list.append(i)
                total_cds += sequence

    gene_list.insert(0, "nonsense|total_cds\n" + total_cds)

    #Obtain UTR list and generate total UTR string, which is useful for later functions
    sequence_list1, sequence_list2, total_sequence = get_sequence_stuff(genes_c)
    seq_lengths1 = [float(len(i)) for i in sequence_list1]
    seq_lengths2 = [float(len(i)) for i in sequence_list2]
    average_length1 = np.mean(seq_lengths1)
    average_length2 = np.mean(seq_lengths2)

    #Get observed and expected frequencies
    observed_counts1 = get_codon_freqs(sequence_list1)
    observed_counts2  = get_codon_freqs(sequence_list2)
    print ('observed frequencies done')

    expected_counts = get_simulations(total_sequence, average_length1, 0)
    print ('expected frequencies done')

    observed_counts1.insert(0, '+1_Frameshift')
    observed_counts2.insert(0, '+2_Frameshift')
    expected_counts.insert(0, 'Simulations')

    csv_total.append(observed_counts1)
    csv_total.append(observed_counts2)
    csv_total.append(expected_counts)

    #Set headers for the CSV
    headers = ["Sequence", 'tag', 'tga', 'gaa', 'gat', 'gta', 'gtt', 'aga', 'agt', 'tgt', 'aag', 'atg', 'ttg', 'n']

    create_csv(headers, csv_total, "tag_test2_downstream.csv")


### FUNCTIONS ###

def get_sequence_stuff(list):

    ''' Obtain the UTR sequence for each gene from the FASTA file '''

    utr_list1 = []
    utr_list2 = []
    total_utr = ''

    #For each gene in the FASTA file...
    for i in list:

        #First generate UTR list, subnested for each codon
        a = i.split("\n")
        utr_seq = "".join(a[1:]).lower()
        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        plus_1 = [utr_seq[i:i+3] for i in range(1, len(utr_seq), 3)]
        plus_2 = [utr_seq[i:i+3] for i in range(2, len(utr_seq), 3)]
        # if codon_seq[0] == 'atg':
        #     if codon_seq[-1] == 'taa' or codon_seq[-1] == 'tga' or codon_seq[-1] == 'tag':
        utr_list1.append(plus_1)
        utr_list2.append(plus_2)

        #Now create total UTR list which will be useful later
        total_utr += utr_seq

    return utr_list1, utr_list2, total_utr


def get_codon_freqs(sequences):

    tag = 0
    tga = 0
    gaa = 0
    gat = 0
    gta = 0
    gtt = 0
    aga = 0
    agt = 0
    tgt = 0
    aag = 0
    atg = 0
    ttg = 0

    n = 0

    for seq in sequences:

        for i in seq:

            codon = i.upper()

            n += 1

            if codon == 'TAG':
                tag += 1
            elif codon == 'TGA':
                tga += 1
            elif codon == 'GAA':
                gaa += 1
            elif codon == 'GAT':
                gat += 1
            elif codon == 'GTA':
                gta += 1
            elif codon == 'GTT':
                gtt += 1
            elif codon == 'AGA':
                aga += 1
            elif codon == 'AGT':
                agt += 1
            elif codon == 'TGT':
                tgt += 1
            elif codon == 'AAG':
                aag += 1
            elif codon == 'ATG':
                atg += 1
            elif codon == 'TTG':
                ttg += 1

    counts = [tag, tga, gaa, gat, gta, gtt, aga, agt, tgt, aag, atg, ttg, n]

    return counts


def get_oof_stops(seq):

    codon_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    oof_codon_seq = [seq[i:i+3] for i in range(1, len(seq), 3)]

    oof_count = 0
    for i in oof_codon_seq:
        if i == 'TAA' or i == 'TGA' or i == 'TAG':
            oof_count += 1

    n = len(oof_codon_seq)
    freq = oof_count / n

    return oof_count, n, freq


def get_simulations(total_sequence, av_length, reading_frame):

    ''' Simulate 10,000 null sequences for each gene based upon dinucleotide content '''

    #Probability of first base - dictionary
    nt_dict = {}
    length = len(total_sequence)

    for i in range(length):
        base = total_sequence[i]
        if base not in nt_dict:
            nt_dict[base] = 1
        else:
            nt_dict[base] += 1

    nt_freq = {k: v / length for k, v in nt_dict.items()}

    print ('nucleotide dictionary done')

    #Second base / Next base
    trans = {}
    for i in range(len(total_sequence)-1):

        dinuc = total_sequence[i:i+2]
        first_dinuc = dinuc[0]
        sec_dinuc = dinuc[1]

        if first_dinuc not in trans:
            trans[first_dinuc] = [sec_dinuc]
        else:
            trans[first_dinuc] += sec_dinuc

    print ('dinucleotide dictionary done')

    #Generate simulations
    total_sequences = []
    total_codons = []
    count = 0

    #This while loop determines how many simulations will be completed
    while count < 1000:

        sim_n = []
        codon_list = []

        sim = []

        np.random.seed()

        #Calculate next base
        first_base = np.random.choice(['a', 'c', 'g', 't'], p=[nt_freq['a'], nt_freq['c'], nt_freq['g'], (1 - nt_freq['a'] - nt_freq['c'] - nt_freq['g'])])
        sim.append(first_base)

        #For one gene simulation
        while (len(sim) < 4):

            #Generate next base
            prev_base = sim[-1]
            next_base = np.random.choice(trans[prev_base], replace=True)
            sim.append(next_base)

        sim = "".join(sim)
        codon_seq = [sim[i:i+3] for i in range(int(reading_frame), len(sim), 3)]
        print (codon_seq)
        total_sequences.append(sim)
        total_codons.append(codon_seq)
        count += 1
        print ('+1 simulation', count, '/10,000')


    freqs = get_codon_freqs(total_codons)

    return freqs


def create_csv(headers, csv_total, filename):

    ''' Output to CSV file '''

    subdir = "TAG_Tests/Test2/CSV"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)


### RUN ###

if __name__ == '__main__':
    main()
