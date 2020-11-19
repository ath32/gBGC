''' Input tsv containing de novo mutations. Convert to mutational matrix table. '''

## IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict
import csv

### CHOOSE SOURCE FOLDER ###

source = './bin10_rec.tsv'
chromosomes = './nature24018-s2/Chromosomes'

### MAIN ###

def main():

    #Set CSV
    csv_total = []
    csv_total.append(['REF', 'ALT', 'Count', 'Total'])

    #Counts
    a_to_c = 0
    a_to_g = 0
    a_to_t = 0

    c_to_a = 0
    c_to_g = 0
    c_to_t = 0

    g_to_a = 0
    g_to_c = 0
    g_to_t = 0

    t_to_a = 0
    t_to_c = 0
    t_to_g = 0

    #Load data
    raw_table = open(os.path.abspath(source)).read()
    tableSplit = [i.split('\t') for i in raw_table.split('\n')[1:] if i != '']
    print (tableSplit[3:])

    #Iterate through each row to get mutational counts
    for row in tableSplit:

        ref = row[2]
        alt = row[3]

        if ref == 'A' and alt == 'C':
            a_to_c += 1
        elif ref == 'A' and alt == 'G':
            a_to_g += 1
        elif ref == 'A' and alt == 'T':
            a_to_t += 1
        elif ref == 'C' and alt == 'A':
            c_to_a += 1
        elif ref == 'C' and alt == 'G':
            c_to_g += 1
        elif ref == 'C' and alt == 'T':
            c_to_t += 1
        elif ref == 'G' and alt == 'A':
            g_to_a += 1
        elif ref == 'G' and alt == 'C':
            g_to_c += 1
        elif ref == 'G' and alt == 'T':
            g_to_t += 1
        elif ref == 'T' and alt == 'A':
            t_to_a += 1
        elif ref == 'T' and alt == 'C':
            t_to_c += 1
        elif ref == 'T' and alt == 'G':
            t_to_g += 1

    #Now get the total counts from the reference genome

    total_a = 0
    total_c = 0
    total_g = 0
    total_t = 0

    # for root, dirs, filenames in os.walk(chromosomes):
    #     for f in filenames:
    #
    #         path = os.path.join(chromosomes, f)
    #         raw_fasta = open(path).read()
    #
    #         genomeSplit = raw_fasta.split('\n')
    #         sequence = "".join(genomeSplit[1:])
    #
    #         total_a += sequence.count('A')
    #         total_c += sequence.count('C')
    #         total_g += sequence.count('G')
    #         total_t += sequence.count('T')

    # Append to CSV
    csv_total.append(['A', 'C', a_to_c, total_a])
    csv_total.append(['A', 'G', a_to_g, total_a])
    csv_total.append(['A', 'T', a_to_t, total_a])
    csv_total.append(['C', 'A', c_to_a, total_c])
    csv_total.append(['C', 'G', c_to_g, total_c])
    csv_total.append(['C', 'T', c_to_t, total_c])
    csv_total.append(['G', 'A', g_to_a, total_g])
    csv_total.append(['G', 'C', g_to_c, total_g])
    csv_total.append(['G', 'T', g_to_t, total_g])
    csv_total.append(['T', 'A', t_to_a, total_t])
    csv_total.append(['T', 'C', t_to_c, total_t])
    csv_total.append(['T', 'G', t_to_g, total_t])

    #Create CSV
    filename = "mutational_counts_bin10_rec.csv"
    subdir = "nature24018-s2/CSV"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        for k in csv_total:
            writer.writerow(k)

### RUN ###

if __name__ == '__main__':
    main()
