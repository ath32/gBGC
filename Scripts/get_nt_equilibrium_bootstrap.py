''' Obtain nucleotide equilibrium frequencies, with bootstrapped error, given a mutational matrix '''

## IMPORTS ###

import os
import re
import numpy as np
import csv
import random
import scipy.stats

### SOURCE GENOME ###

source = 'nature24018-s2/CSV/mutational_counts_bin10_rec.csv'

### MAIN ###

def main():

    #First, unpack and get switch frequencies
    raw_file = open(source).read()
    split = raw_file.split('\n')
    counts = [float(i.split(',')[2]) for i in split[1:] if i != '']
    totals = [float(i.split(',')[3]) for i in split[1:] if i != '']
    freqs = [x/y for x, y in zip(counts, totals)]

    #Observed freqs
    a_eq, t_eq, g_eq, c_eq, gc_eq = get_eq(freqs)

    print (a_eq, c_eq, g_eq, t_eq, gc_eq*100)

    #Bootstrap!

    #First, get a list of all mutations
    mutations = ['a_c', 'a_g', 'a_t', 'c_a', 'c_g', 'c_t', 'g_a', 'g_c', 'g_t', 't_a', 't_c', 't_g']
    sampling_list = []
    for i,n in enumerate(mutations):
        name = n
        count = counts[i]
        no = 0
        while no < count:
            sampling_list.append(name)
            no += 1

    #Next, sample with replacement a new list
    bootstrapped_gc = []
    counter = 0
    while counter < 100:
        sample = random.choices(sampling_list, k=len(sampling_list))
        new_counts = [sample.count('a_c'), sample.count('a_g'), sample.count('a_t'),
                        sample.count('c_a'), sample.count('c_g'), sample.count('c_t'),
                        sample.count('g_a'), sample.count('g_c'), sample.count('g_t'),
                        sample.count('t_a'), sample.count('t_c'), sample.count('t_g')]
        frequencies = [x/y for x, y in zip(new_counts, totals)]
        a_eq, t_eq, g_ew, c_eq, gc_eq = get_eq(frequencies)
        bootstrapped_gc.append(gc_eq)
        counter += 1

    mean, lower, upper = mean_confidence_interval(bootstrapped_gc, confidence=0.95)
    print (mean, lower, upper)

def get_eq(freqs):

    a_c = freqs[0]
    a_g = freqs[1]
    a_t = freqs[2]
    c_a = freqs[3]
    c_g = freqs[4]
    c_t = freqs[5]
    g_a = freqs[6]
    g_c = freqs[7]
    g_t = freqs[8]
    t_a = freqs[9]
    t_c = freqs[10]
    t_g = freqs[11]

    #Make equations
    a_coef1 = 1 - (1 - a_t - a_g - a_c) + c_a
    t_coef1 = c_a - t_a
    g_coef1 = c_a - g_a
    constant1 = c_a

    t_coef2 = 1 - (1 - t_a - t_g - t_c) + c_t
    a_coef2 = c_t - a_t
    g_coef2 = c_t - g_t
    constant2 = c_t

    g_coef3 = 1 - (1 - g_a - g_t - g_c) + c_g
    a_coef3 = c_g - a_g
    t_coef3 = c_g - t_g
    constant3 = c_g

    #Solve
    equations = np.array([[a_coef1, t_coef1, g_coef1], [a_coef2, t_coef2, g_coef2], [a_coef3, t_coef3, g_coef3]])
    answer = np.array([constant1, constant2, constant3])

    #Solve the simultaneous equations
    solved = np.linalg.solve(equations,answer)
    a_eq = solved[0]
    t_eq = solved[1]
    g_eq = solved[2]
    c_eq = 1 - a_eq - t_eq - g_eq

    #Print results
    gc_eq = g_eq + c_eq

    return a_eq, t_eq, g_eq, c_eq, gc_eq

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

if __name__ == '__main__':
  main()
