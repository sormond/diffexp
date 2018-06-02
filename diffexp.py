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

""" def error_calling(df) :"""

def merge_df(file1, file2) :
    f1 = loadf(file1)
    f2 = loadf(file2)
    #error_calling(f1)
    #error_calling(f2)
    x = []
    x.append(f1.iloc[1])
    df = merge(f1, f2)
    """ for i in range(0, len(f1)) :
        x = (f1.iloc[i])
        if type(x) is not int or type(x) is not float :
            print("Warning: In row number %s of dataframe, expression value is invalid. Please ensure a number value is present" % i) """
    return df

import numpy as np
import scipy
from scipy.stats import chi2_contingency

# this function works, but need to find out how to get it to return 'nan' when two exp. values are '0'
def chisqtest(df) :
    A = (sum(merged['sampleA']))
    B = (sum(merged['sampleB']))
    xsquared_values = []
    p_values =[]
    np_df = df.as_matrix()
    for i in range(0, len(np_df)) :
        a = np_df[i,0]
        b = np_df[i,1]
        # this is a quick fix for 0 values in the array, but need to find workaround for this
        if a == 0 and b == 0 :
            xsquared_values.append('nan')
            p_values.append('nan')
        else :   
            x = np.array([[a,b], [A, B]])
            chi = chi2_contingency(x)
            xsquared_values.append(chi[0])
            p_values.append(chi[1])

#includes both'merge_df' and 'chisqtest' function in one
def stats_dataframe(file1, file2) :
    f1 = loadf(file1)
    f2 = loadf(file2)
    #error_calling(f1)
    #error_calling(f2)
    mdf = merge(f1, f2)
    A = (sum(mdf['sampleA']))
    B = (sum(mdf['sampleB']))
    xsquared_values = []
    p_values =[]
    np_df = mdf.as_matrix()
    for i in range(0, len(np_df)) :
        a = np_df[i,0]
        b = np_df[i,1]
        # don't know if it's correct to have this as 'NaN'?
        if a == 0 and b == 0 :
            xsquared_values.append('nan')
            p_values.append('nan')
        else :   
            x = np.array([[a,b], [A, B]])
            chi = chi2_contingency(x)
            xsquared_values.append(chi[0])
            p_values.append(chi[1])
    xsquared_values = pd.Series(xsquared_values)
    p_values = pd.Series(p_values)
    mdf['x_sq_values'] = xsquared_values.values
    mdf['p_values'] = p_values.values
    return mdf
        


import numpy as np
    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--sample1', help = "File of sample 1.")
parser.add_argument('-2', '--sample2', help = "File of sample 2.")
""" parser.add_argument('-o', '--outstem', help = "Stem name for output file. Default stem is 'output' if not given.") """

args = parser.parse_args()

""" merged = merge_df(args.sample1, args.sample2)
print(merged) """
""" print(merged.iloc[100]) """

import scipy
from scipy.stats import chi2_contingency
""" obs = np.array(np_df).T
print(obs.shape)
print(chisquare(obs)) """
""" print(sum(merged['sampleA']))
print(sum(merged['sampleB'])) """
""" obs = np.array([[63, 14], [15300, 5797]])
print(chi2_contingency(obs))
print(merged.iloc[0,1]) """
""" A = 100
B = 200 """
""" x = np.array([(pd.to_numeric(merged.iloc[5,0])), (pd.to_numeric(merged.iloc[4,1]))], [(100), (200)])"""
""" x = ((np_df[1]), (np_df[2]))
y = (np_df[1,0])
print(y) """
dataframe = stats_dataframe(args.sample1, args.sample2)
print(dataframe)

#Below code works
""" p = 1
q = np.array([[p, p], [p, p]])
print(chi2_contingency(q))
print(q) """



""" df = pd.DataFrame(np.random.rand(50, 4), columns=['a', 'b', 'c', 'd'])
df.plot.scatter(x='a', y='b')
plt.show()
import matplotlib.pyplot as plt
plt.scatter(sample1[0], sample1[1])
plt.show() """
""" To plot multiple column groups in a single axes, repeat plot method specifying target ax. It is recommended to specify color and label keywords to distinguish each groups.
In [62]: ax = df.plot.scatter(x='a', y='b', color='DarkBlue', label='Group 1');
In [63]: df.plot.scatter(x='c', y='d', color='DarkGreen', label='Group 2', ax=ax);
"""