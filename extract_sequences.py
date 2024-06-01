# Author: Henrik Seidel
# Date: 2018-07-03
# This script extracts all sequences from a fasta file into single sequences

import sys
import os.path
import gzip
import re

if len(sys.argv) != 3:
    sys.exit("ERROR: you have to specify the name of the fasta file " +
    "and the name of the output directory")

file_name = sys.argv[1]
out_dir = sys.argv[2]

# Assert we have a filename as first argument
if not(os.path.isfile(file_name)):
    sys.exit("ERROR: No such file: \"" + file_name + "\"")

# Assert the name of the output directory is okay
if not(os.path.exists(out_dir)):
    os.mkdir(out_dir)
elif not(os.path.isdir(out_dir)):
    sys.exit("ERROR: \"" + out_dir + "\" exists but is not a directory")

# Choose function for opening file, depending on compression status of file
if (file_name.endswith(".gz")):
    this_open = gzip.open
    file = gzip.open(file_name, mode = "rt")
else:
    this_open = open

n = 0
with this_open(file_name, mode = "rt") as file:
    out_file = None
    try:
        for line in file:
            if line.startswith(">"):
                if out_file is not None:
                    out_file.close()
                seq_name = re.match("^>([^ ]+)", line)[1]
                n += 1
                out_file_name = '{:s}/{:03d}_{:s}.fa'.format(out_dir, n, seq_name)
                print('Writing {:s} to {:s}'.format(seq_name, out_file_name))
                out_file = open(out_file_name, mode = "wt")
            out_file.write(line)
    finally:
        out_file.close()
