#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v0.1
# DATE        :  2021-4-30
# DESCRIPTION :  Calculate jaccard distance between all CRISPR arrays and lookup core-genome SNP differences between isolates encoding those arrays. Plot as scatterplot.

import sys
import argparse
import textwrap as _textwrap


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
