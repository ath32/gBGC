''' Script to get stop codon equilibrium frequencies, with bootstrapped error,
    given a mutational matrix '''

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
    taa_eq, tga_eq, tag_eq = get_eq(freqs)
    print (taa_eq, tga_eq*100, tag_eq)

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
    bootstrapped_taa = []
    counter = 0
    while counter < 100:
        sample = random.choices(sampling_list, k=len(sampling_list))
        new_counts = [sample.count('a_c'), sample.count('a_g'), sample.count('a_t'),
                        sample.count('c_a'), sample.count('c_g'), sample.count('c_t'),
                        sample.count('g_a'), sample.count('g_c'), sample.count('g_t'),
                        sample.count('t_a'), sample.count('t_c'), sample.count('t_g')]
        frequencies = [x/y for x, y in zip(new_counts, totals)]
        taa_eq, tga_eq, tag_eq = get_eq(frequencies)
        bootstrapped_taa.append(tag_eq)
        counter += 1

    mean, lower, upper = mean_confidence_interval(bootstrapped_taa, confidence=0.95)
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
    taa_coef1 = 1 - (1 - a_g - a_g) + g_a
    tga_coef1 = g_a - g_a
    constant1 = g_a

    tga_coef2 = 1 - (1 - g_a - (2*g_a*a_g)) + (2*g_a*a_g)
    taa_coef2 = (2*g_a*a_g) - a_g
    constant2 = (2*g_a*a_g)

    equations = np.array([[taa_coef1, tga_coef1], [taa_coef2, tga_coef2]])
    answer = np.array([constant1, constant2])

    #Solve the simultaneous equations
    solved = np.linalg.solve(equations,answer)
    taa_eq = solved[0]
    tga_eq = solved[1]
    tag_eq = 1 - taa_eq - tga_eq

    #Print results
    # print (taa_eq, tga_eq, tag_eq)
    # #Check the results
    # print (np.allclose(np.dot(equations, solved), answer))

    return taa_eq, tga_eq, tag_eq

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

if __name__ == '__main__':
  main()
