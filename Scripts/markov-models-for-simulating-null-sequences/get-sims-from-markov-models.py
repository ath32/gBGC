### IMPORT STATEMENTS ###

import os
import numpy as np
import csv
import scipy.stats as st
import random

### SOURCE FILES ###

dinuc_source = './dinuceotide_equilibria.csv'

### MAIN ###

def main():

    #Empty lists to contain results
    results = []
    gc_results = []

    #For 100 sets of sims for example, increase this as needed
    for i in range(0,100):

        #Load equilibria to be entered into the Markov models
        dinuc_eq = eq_from_csv(dinuc_source)

        #Simulate sequences - input equilibria and number of sims to be simulated
        dinuc_sims = markov_sims(dinuc_eq, 100)

        #Calculate trinucleotide frequencies from simulation
        dinuc_trinucs = get_trinucleotides(dinuc_sims)

        #Calculate GC content from simulation
        gc = get_gc(dinuc_sims)

        #Calculate TAA, TGA and TAG equilibrium frequencies
        taa = [i[1] for i in dinuc_trinucs if i[0] == 'TAA']
        tga = [i[1] for i in dinuc_trinucs if i[0] == 'TGA']
        tag = [i[1] for i in dinuc_trinucs if i[0] == 'TAG']
        rel_tga = tga[0] / (taa[0] + tga[0] + tag[0])

        #Append results
        results.append(rel_tga)
        gc_results.append(gc)

    #Print results - mean and std of simulations
    print ('TGA mean:', np.mean(results))
    print ('TGA standard dev:', np.std(results))

    print ('GC mean:', np.mean(gc_results))
    print ('GC standard dev:', np.std(gc_results))

### FUNCTIONS ###

def get_gc(sims):

    gc_list = []
    for i in sims:
        gc = (i.count('G') + i.count('C')) / len(i)
        gc_list.append(gc)

    return np.mean(gc_list)

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

    #Convert the frequencies into dictionaries for selecting base according to prev base
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

    #Count is the number of simulations to be generated
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

        #Join into string
        sim = "".join(sim)

        #Append simulated sequence to list of sims
        total_sequences.append(sim)
        count += 1
        print ('+1 simulation', count, '/1000')

    return total_sequences

def get_trinucleotides(seqs):

    dict = {}
    n = 0

    #Extract all trinucleotides in all reading frames
    for sequence in seqs:
        trinucleotides = []
        for i in range(len(sequence)-2):
            trinucleotide = sequence[i:i+3]
            if 'N' not in trinucleotide:
                trinucleotides.append(trinucleotide)

        #Then count the number of each type
        for trinuc in trinucleotides:
            n += 1
            if trinuc in dict:
                dict[trinuc] += 1
            else:
                dict[trinuc] = 1

    #Calculate frequencies
    freq_dict = {k: v / n for k, v in dict.items()}
    dictlist = []
    for key, value in freq_dict.items():
        entry = [key,value]
        dictlist.append(entry)

    return dictlist

### RUN ###

if __name__ == '__main__':
    main()
