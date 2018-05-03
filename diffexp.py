#!/usr/bin/env python3

""" wkdirectory = "/Users/sormond/Desktop//diffexp/DataFiles" # add working directory here
import os
os.chdir(wkdirectory) """
import pandas as pd
def loadf(file) :
    df = pd.read_table(file, index_col = 0)
    return df
    """ with open(file, 'r') as r:
        df = r.readlines()
        return list """

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--sample1', help = "File of sample 1.")
parser.add_argument('-2', '--sample2', help = "File of sample 2.")
""" parser.add_argument('-o', '--outstem', help = "Stem name for output file. Default stem is 'output' if not given.") """

args = parser.parse_args()

sample1 = loadf(args.sample1)
sample2 = loadf(args.sample2)

import matplotlib.pyplot as plt
plt.scatter(sample1[0], sample1[1])
plt.show()
