#!/usr/bin/env python3

# Make a dataframe from .txt files. Index = gene ID. One column = expression value

""" import matplotlib.pyplot as plt """
import pandas as pd
""" import numpy as np """

""" wkdirectory = "/Users/sormond/Desktop//diffexp/DataFiles" # add working directory here
import os
os.chdir(wkdirectory) """

pd.options.display.max_rows = 999

def loadf(file) :
    df = pd.read_csv(file, sep = "\t")
    l = df['gene'].tolist()
    df.index = l
    del df['gene']
    return df

def merge(df1, df2) :
    df = df1.join(df2, how='outer')
    return df

def merged_df(file1, file2) :
    f1 = loadf(file1)
    f2 = loadf(file2)
    df = merge(f1, f2)
    return df

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--sample1', help = "File of sample 1.")
parser.add_argument('-2', '--sample2', help = "File of sample 2.")
""" parser.add_argument('-o', '--outstem', help = "Stem name for output file. Default stem is 'output' if not given.") """

args = parser.parse_args()

merged = merged_df(args.sample1, args.sample2)
print(merged)

""" df = pd.DataFrame(np.random.rand(50, 4), columns=['a', 'b', 'c', 'd'])
df.plot.scatter(x='a', y='b')
plt.show() """
""" import matplotlib.pyplot as plt
plt.scatter(sample1[0], sample1[1])
plt.show() """

""" To plot multiple column groups in a single axes, repeat plot method specifying target ax. It is recommended to specify color and label keywords to distinguish each groups.
In [62]: ax = df.plot.scatter(x='a', y='b', color='DarkBlue', label='Group 1');
In [63]: df.plot.scatter(x='c', y='d', color='DarkGreen', label='Group 2', ax=ax);
"""