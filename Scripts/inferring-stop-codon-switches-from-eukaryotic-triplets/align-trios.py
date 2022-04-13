### IMPORT STATEMENTS ###

import os
from generic import run_in_parallel

### SOURCE FILES ###

source_in = "./Primates/Unaligned"
source_out = "./Primates/Aligned"

### MAIN ###

def main():

    filenames = get_files(source_in)
    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(filenames, ['foo', source_in], command, workers=workers)

    csv_total = []

    for process in processes:
        output = process.get()

### FUNCTIONS ###

def get_files(source):

    files = []
    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            files.append(f)

    return files

def command(filenames, source_in):

    for root, dirs, filenames in os.walk(source_in):
        for f in filenames:
            #Define path for input/output files
            input = os.path.join(source_in, f)
            output = os.path.join(source_out, f)
            #Send input file to MAFFT, recieve aligned file
            cmd = "mafft --auto %s > %s" % (input, output)
            os.system(cmd)

### RUN ###

if __name__ == '__main__':
    main()
