### IMPORT STATEMENTS ###

import os
import generic as gen
import re
import collections
import copy
import numpy as np
import random
import shutil
from pathlib import Path
from tqdm import tqdm

### SOURCE FILES - Ingroup 1 FASTA, Ingroup 2 FASTA, Outgroup FASTA, Orthologs list as TSV ###

in1_source = './FASTA/mus_musculus.fa'
in2_source = './FASTA/mus_spretus.fa'
out1_source = './FASTA/rattus_norvegicus.fa'
ortho_source = './Ortho/mouse_orthologs.txt'

### MAIN ###

def main():

    #Load orthologs
    raw_ortho = open(ortho_source).read()
    ortho_split = raw_ortho.split('\n')

    #Extract all orthologs
    for gene in tqdm(ortho_split[1:]):
        #Get ids for each genome
        sections = gene.split('\t')
        #Ignore all non-1:1 orthologs
        if sections[0] != '' and sections[1] != '' and sections[2] != '':

            in1_id = sections[0]
            in2_id = sections[1]
            out_id = sections[2]

            file = ''

            #Extract the relevant in1 gene
            raw_in1 = open(in1_source).read()
            in1_split = raw_in1.split('>')
            for a in in1_split:
                if a != '':
                    a_split = a.split('\n')
                    line = a_split[0]
                    line_split = line.split('|')
                    gene = line_split[0].replace('ID=gene:', '')
                    sequence = "".join(a_split[1:])
                    if in1_id in line and sequence != 'Sequence unavailable':
                        if in1_id not in file:
                            file += '>' + line + '\n' + sequence + '\n'

            #Extract the relevant in2 gene
            raw_in2 = open(in2_source).read()
            in2_split = raw_in2.split('>')
            for b in in2_split:
                if b != '':
                    b_split = b.split('\n')
                    line = b_split[0]
                    line_split = line.split('|')
                    gene = line_split[0].replace('ID=gene:', '')
                    sequence = "".join(b_split[1:])
                    if in2_id in line and sequence != 'Sequence unavailable':
                        if in2_id not in file:
                            file += '>' + line + '\n' + sequence + '\n'

            #Extract the relevant out gene
            raw_out = open(out1_source).read()
            out_split = raw_out.split('>')
            for c in out_split:
                if c != '':
                    c_split = c.split('\n')
                    line = c_split[0]
                    line_split = line.split('|')
                    gene = line_split[0].replace('ID=gene:', '')
                    sequence = "".join(c_split[1:])
                    if out_id in line and sequence != 'Sequence unavailable':
                        if out_id not in file:
                            file += '>' + line + '\n' + sequence + '\n'

            #Write output file
            if len(file.split('>')) > 3:
                out_source = "./Unaligned"
                out_path = os.path.join(out_source, in1_id + '.fa')
                if os.path.isfile(out_path) == False:
                    out_file = open(out_path, 'a')
                    out_file.write(file)
                    out_file.close()

### RUN ###

if __name__ == '__main__':
    main()
