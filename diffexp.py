#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--sample1', help = "File of sample 1.")
parser.add_argument('-2', '--sample2', help = "File of sample 2.")
parser.add_argument('-o', '--outstem', help = "Stem name for output file. Default stem is 'output' if not given.")

args = parser.parse_args()

sample_1 = args.sample1
sample_2 = args.sample2

def load(f)
    with open(f, 'r') as f:
        fastq_list = f.readlines()

