### IMPORT STATEMENTS ###

import os
import re
import numpy as np
import csv
from shutil import copyfile

### SOURCE FILES - Table of IDs with GC contents/Recombination rates, Alignments folder, Directories for GC-rich/GC-poor or HRGs/LRGs  ###

table_source = './Primates/CSV/recombination_ids.csv'
alignment_source = './Primates/Aligned'
source1 = './Primates/Aligned_HRGs'
source2 = './Primates/Aligned_LRGs'

### MAIN ###

def main():

    #Get list of IDs for which we have an alignment
    ids = []
    for root, dirs, filenames in os.walk(alignment_source):
        for f in filenames:
            path = os.path.join(alignment_source, f)
            raw_file = open(path).read()
            file_ids = []
            raw_split = raw_file.split('>')
            for i in raw_split:
                if i != '':
                    splt = i.split('\n')
                    id = splt[0].split(';')[0].replace('ID=gene:', '')
                    file_ids.append(id)
            ids.append(file_ids[1])

    #Empty lists to hold GC-rich/GC-poor or HRG/LRG gene ids
    one = []
    two = []

    #Empty lists to hold the GC contents/recombination rates in each gene bin
    oneNumbers = []
    twoNumbers = []

    #Find the 50% cut off for gene binning by GC content or recombination rate
    table = open(os.path.abspath(table_source)).read()
    rows = table.split('\n')
    values_list = [i.split(',')[:2] for i in rows[1:] if i != '']
    values = [float(i[1]) for i in values_list if i[0] in ids and float(i[1]) > 0]
    q1 = np.quantile(values, 0.50)

    #Assign aligned orthologs to gene bins
    valid = 0
    for row in rows[1:]:
        if row != '':
            row_split = row.split(',')
            if row_split[0] in ids:
                valid += 1
                if 0 < float(row_split[1]) <= q1:
                    one.append(row_split[0])
                    oneNumbers.append(float(row_split[1]))
                elif q1 < float(row_split[1]):
                    two.append(row_split[0])
                    twoNumbers.append(float(row_split[1]))

    #Calculate mean GC content/recombination rate of each bin
    print (np.mean(oneNumbers))
    print (np.mean(twoNumbers))

    #Filter aligned files and move them to their new directory
    for root, dirs, filenames in os.walk(alignment_source):
        for f in filenames:

            #Grab gene ids in the aligned file
            path = os.path.join(alignment_source, f)
            raw_file = open(path).read()
            file_ids = []
            raw_split = raw_file.split('>')
            for i in raw_split:
                if i != '':
                    splt = i.split('\n')
                    id = splt[0].split(';')[0].replace('ID=gene:', '')
                    file_ids.append(id)

            #Copy file to new directory
            for gene in file_ids:
                if gene in one:
                    destination = os.path.join(source1, f)
                    copyfile(path, destination)
                elif gene in two:
                    destination = os.path.join(source2, f)
                    copyfile(path, destination)

### RUN ###

if __name__ == '__main__':
    main()
