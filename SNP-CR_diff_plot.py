#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v0.1
# DATE        :  2021-4-30
# DESCRIPTION :  Calculate jaccard similarity between all CRISPR arrays and lookup core-genome SNP differences between isolates encoding those arrays. Plot as scatterplot.

import sys
import argparse
import pickle
import gzip
from itertools import combinations
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="Calculate jaccard similarity between all CRISPR arrays and lookup core-genome SNP differences between isolates encoding those arrays. Plot as scatterplot.")
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

Jaccard_list = []
Core_SNP_list = []

# Add the SNP distances of isolates encoding each array first (i.e. how far apart are isolates with the same array; jaccard similarity = 1)

for array, reps in array_rep_dict.items():
    reps = list(set(reps)) # Use set in case any of the arrays are duplicated in a genome
    if len(reps) > 1:
        for combo in combinations(reps, 2):
            if combo in SNP_diff_dict.keys():
                Core_SNP_list.append(SNP_diff_dict[combo])
            elif tuple(reversed(combo)) in SNP_diff_dict.keys():
                Core_SNP_list.append(SNP_diff_dict[tuple(reversed(combo))])
            else:
                print("Can't find {}".format(combo))
                continue
            Jaccard_list.append(1)

plt.hist(Core_SNP_list, density=False, bins=1000)
plt.title("Histogram of distribution of SNP distances\nbetween isolates encoding identical arrays")
plt.yscale("log")
plt.xlabel('Number of SNPs between isolates encoding an identical array')
plt.ylabel('Bin count (log10)')
plt.tight_layout()
plt.savefig(args.out_prefix + 'identical_array_snp_distances.png', dpi=300)
plt.close()
