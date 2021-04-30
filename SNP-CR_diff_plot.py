#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v0.1
# DATE        :  2021-4-30
# DESCRIPTION :  Calculate jaccard distance between all CRISPR arrays and lookup core-genome SNP differences between isolates encoding those arrays. Plot as scatterplot.

import sys
import argparse
import pickle
import gzip
from itertools import combinations

parser = argparse.ArgumentParser(
    description="Calculate jaccard distance between all CRISPR arrays and lookup core-genome SNP differences between isolates encoding those arrays. Plot as scatterplot.")
parser.add_argument(
    "-a", dest="array_representatives", required = True,
    help="Array_ID_representatives file output by minced2network.py."
    )
parser.add_argument(
    "-n", dest="array_network", required = True,
    help="array network file output by minced2network.py."
    )
parser.add_argument(
    "-d", dest="diff_dict", required = True,
    help="Pickled SNP diff dict output by snps2diffdict.py."
    )

parser.add_argument(
    "-z", action='store_true',
        help="Indicate that the diff dict is gzipped."
    )
parser.add_argument(
    "-o", dest="out_prefix", required = True,
    help="Path and file prefix for output file(s)."
    )

args = parser.parse_args(sys.argv[1:])

# SNP_diff_dict of the format {('isolate1', 'isolate2') : num_SNPs}

if args.z:
    with gzip.open(args.diff_dict, 'rb') as fin:
        SNP_diff_dict = pickle.load(fin)
else:
    with open(args.diff_dict, 'rb') as fin:
        SNP_diff_dict = pickle.load(fin)

array_rep_dict = {} # {'Array_ID' : ['Isolate1', 'Isolate2']} where the list includes all isolates with that array

with open(args.array_representatives, 'r') as fin:
    for line in fin.readlines()[1:]:
        elements = line.split()
        array_rep_dict[elements[0]] = elements[1:]
