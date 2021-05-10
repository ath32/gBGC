### IMPORTS ###

import os
import numpy as np
import csv
import scipy.stats as st
import random

### SOURCE ###

dinuc_source = './Equilibrium/Di/Bespoke_totals/New/dinuc_eq_9.csv'

### MAIN ###

def main():

    #Get equilibria to be entered into the Markov models
    dinuc_eq = eq_from_csv(dinuc_source)
    dinuc_sims = markov_sims(dinuc_eq, 1000)
    dinuc_trinucs = get_trinucleotides(dinuc_sims, 'Dinuc')
    taa = [i[1] for i in dinuc_trinucs if i[0] == 'TAA']
    tga = [i[1] for i in dinuc_trinucs if i[0] == 'TGA']
    tag = [i[1] for i in dinuc_trinucs if i[0] == 'TAG']
    rel_tga = tga[0] / (taa[0] + tga[0] + tag[0])
    print (rel_tga)


### FUNCTIONS ###

def eq_from_csv(source):

    raw_file = open(os.path.abspath(source)).read()
    rawSplit = raw_file.split('\n')[1:]
    equilibria = []
    for i in rawSplit:
        if i != '':
            line_split = i.split(',')
            dinucleotide = line_split[0]
            freq = line_split[1]
            equilibria.append([dinucleotide, freq])

    return equilibria

def markov_sims(eq, av_length):

    #First, convert the frequencies into dictionaries for selecting base according to prev base
    a_trans = {}
    total_a = 0
    for section in eq:
        dinuc = section[0]
        first_nt = dinuc[0]
        second_nt = dinuc[1]
        freq = float(section[1])
        if first_nt == 'A':
            total_a += freq
            if second_nt not in a_trans:
                a_trans[second_nt] = freq
    for key, value in a_trans.items():
        a_trans[key] = value / total_a

    c_trans = {}
    total_c = 0
    for section in eq:
        dinuc = section[0]
        first_nt = dinuc[0]
        second_nt = dinuc[1]
        freq = float(section[1])
        if first_nt == 'C':
            total_c += freq
            if second_nt not in c_trans:
                c_trans[second_nt] = freq
    for key, value in c_trans.items():
        c_trans[key] = value / total_c

    g_trans = {}
    total_g = 0
    for section in eq:
        dinuc = section[0]
        first_nt = dinuc[0]
        second_nt = dinuc[1]
        freq = float(section[1])
        if first_nt == 'G':
            total_g += freq
            if second_nt not in g_trans:
                g_trans[second_nt] = freq
    for key, value in g_trans.items():
        g_trans[key] = value / total_g

    t_trans = {}
    total_t = 0
    for section in eq:
        dinuc = section[0]
        first_nt = dinuc[0]
        second_nt = dinuc[1]
        freq = float(section[1])
        if first_nt == 'T':
            total_t += freq
            if second_nt not in t_trans:
                t_trans[second_nt] = freq
    for key, value in t_trans.items():
        t_trans[key] = value / total_t

    #Now run sims
    total_sequences = []
    count = 0

    while count < 1000:
        np.random.seed()
        sim = []

        #Start with a random nucleotide
        first_dinuc = np.random.choice([i[0] for i in eq], p=[i[1] for i in eq])
        first_base = first_dinuc[0]
        sim.append(first_base)

        #For one gene simulation
        while (len(sim) <= av_length):

            #Generate next base
            prev_base = sim[-1]
            if prev_base == 'A':
                next_base = np.random.choice(['A', 'C', 'G', 'T'], p = [a_trans['A'], a_trans['C'], a_trans['G'], a_trans['T']])
                sim.append(next_base)
            elif prev_base == 'C':
                next_base = np.random.choice(['A', 'C', 'G', 'T'], p = [c_trans['A'], c_trans['C'], c_trans['G'], c_trans['T']])
                sim.append(next_base)
            elif prev_base == 'G':
                next_base = np.random.choice(['A', 'C', 'G', 'T'], p = [g_trans['A'], g_trans['C'], g_trans['G'], g_trans['T']])
                sim.append(next_base)
            elif prev_base == 'T':
                next_base = np.random.choice(['A', 'C', 'G', 'T'], p = [t_trans['A'], t_trans['C'], t_trans['G'], t_trans['T']])
                sim.append(next_base)

        sim = "".join(sim)
        total_sequences.append(sim)
        count += 1
        print ('+1 simulation', count, '/10,000')

    return total_sequences

def get_trinucleotides(seqs, flag):

    dict = {}
    n = 0

    for sequence in seqs:

        if flag == 'CDS':
            trinucleotides = []
            #Ignore the "stop" codon
            for i in range(len(sequence)-5):
                trinucleotide = sequence[i:i+3]
                if 'N' not in trinucleotide:
                    trinucleotides.append(trinucleotide)
            del trinucleotides[0::3]
            for trinuc in trinucleotides:
                n += 1
                if trinuc in dict:
                    dict[trinuc] += 1
                else:
                    dict[trinuc] = 1

        else:
            trinucleotides = []
            for i in range(len(sequence)-2):
                trinucleotide = sequence[i:i+3]
                if 'N' not in trinucleotide:
                    trinucleotides.append(trinucleotide)
            for trinuc in trinucleotides:
                n += 1
                if trinuc in dict:
                    dict[trinuc] += 1
                else:
                    dict[trinuc] = 1

    freq_dict = {k: v / n for k, v in dict.items()}
    dictlist = []
    for key, value in freq_dict.items():
        entry = [key,value]
        dictlist.append(entry)

    return dictlist

### RUN ###

if __name__ == '__main__':
    main()
