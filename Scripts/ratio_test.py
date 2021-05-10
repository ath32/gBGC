''' Change source file and adjust list indexes to test expression or recombination data '''

### IMPORTS ###

import os
import re
import numpy as np
import csv
import random

### SOURCE ###

source = './Expression/gene_information_bio_switches.csv'

### MAIN ###

def main():

    csv_total = []

    #First, read in the CSV and split into rows
    rawCSV = open(os.path.abspath(source)).read()
    csvSplit = [i.split(',') for i in rawCSV.split('\n')[1:] if i != '']
    gc_contents = [float(i[3]) for i in csvSplit if i != '']
    median = np.median(gc_contents)
    csvSplit = [i for i in csvSplit if float(i[3]) < median]

    #Get upper and lower boundaries for defining HEGs and LEGs
    expression_list = []
    for i in csvSplit:
        if len(i) > 1:
            if float(i[3]) < median:
                expression_list.append(float(i[2]))

    thresholds = [0.5, 0.6, 0.7, 0.8, 0.9]

    count = 0
    for i in thresholds:
        lower_quant = np.quantile(expression_list, i-count)
        higher_quant = np.quantile(expression_list, i)
        all_ids = [i[1] for i in csvSplit if len(i) > 3]
        hegs = [i[1] for i in csvSplit if len(i) > 3 and float(i[2]) >= higher_quant]
        legs = [i[1] for i in csvSplit if len(i) > 3 and float(i[2]) <= lower_quant]

        #Calculate real scores
        a, b, c, d = analyse_switches(hegs, legs, csvSplit)
        conservation_score = (a / b) / (c / d)
        print (i, len(hegs), len(legs), conservation_score, (a/b), (c/d))
        count += 0.2

### FUNCTIONS ###

def analyse_switches(hegs, legs, csvSplit):

    #Calculate a, b, c, and d
    a1 = 0
    a2 = 0
    b1 = 0
    b2 = 0
    c1 = 0
    c2 = 0
    d1 = 0
    d2 = 0

    for i in csvSplit:
        if len(i) > 3:
            id = i[1]
            PxAbundance = float(i[2])
            switch = i[7]
            #HEGs
            if id in hegs:
                if switch == 'TAA_TGA':
                    a1 += 1
                if switch == 'TAA_TGA' or switch == 'TAA_TAG' or switch == 'TAA_stay':
                    b1 += 1
                if switch == 'TGA_TAA':
                    c1 += 1
                if switch == 'TGA_TAA' or switch == 'TGA_TAG' or switch == 'TGA_stay':
                    d1 += 1
            #LEGs
            if id in legs:
                if switch == 'TAA_TGA':
                    a2 += 1
                if switch == 'TAA_TGA' or switch == 'TAA_TAG' or switch == 'TAA_stay':
                    b2 += 1
                if switch == 'TGA_TAA':
                    c2 += 1
                if switch == 'TGA_TAA' or switch == 'TGA_TAG' or switch == 'TGA_stay':
                    d2 += 1

    a = (a1 / b1)
    b = (c1 / d1)
    c = (a2 / b2)
    d = (c2 / d2)

    return a, b, c, d

### RUN ###

if __name__ == '__main__':
    main()
