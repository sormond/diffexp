#!/usr/bin/env python3

""" import matplotlib.pyplot as plt """
import pandas as pd
import numpy as np
import numpy.ma as ma # for masking functions
import scipy
from scipy.stats import chi2_contingency
import statsmodels.stats.multitest
""" from statsmodels.stats import stats
import stats.multitest """

pd.options.display.max_rows = 999 #delete when finalised

def loadf(file) :
    df = pd.read_csv(file, sep = "\t")
    global l
    l = df['gene'].tolist()
    df.index = l
    del df['gene']
    return df

def merge(df1, df2) :
    df = df1.join(df2, how='outer')
    return df

# includes both'merge_df' and 'chisqtest' function in one
# function which creates pandas data frame containing list of genes and expression values
# and adds columns containing calculated chi x-squared values and p-values
def stats_dataframe(file1, file2) :
    # loads in two expression value files and creates dataframes by calling function 'loadf'
    f1 = loadf(file1) 
    f2 = loadf(file2)
    # merges two dataframes into one
    mdf = merge(f1, f2)
    # sums the total of each expession column for statistical calculations
    A = (sum(mdf['sampleA']))
    B = (sum(mdf['sampleB']))
    # creates empty lists for statistical values to be appended to
    xsquared_values = []
    p_values =[]
    # creates numpy matrix from merged dataframe
    np_df = mdf.as_matrix()
    # iterates over each row of matrix to calculate statistical values for each row
    for i in range(0, len(np_df)) :
        a = np_df[i,0]
        b = np_df[i,1]
        """ don't know if it's correct to have this as 'NaN'? """
        # 'if' 'else' statement which prevents error if two expression values for same gene = zero
        if a == 0 and b == 0 :
            xsquared_values.append(np.nan)
            p_values.append(np.nan)
        else :
            # calculates x-squared and p-values using chi contigency table and appends values to previously defined lists
            x = np.array([[a,b], [A, B]])
            chi = chi2_contingency(x)
            xsquared_values.append(chi[0])
            p_values.append(chi[1])
    p_values = np.array(p_values)
    mask = np.isfinite(p_values)
    pval_corrected = np.empty(p_values.shape)
    null_hypothesis = np.empty(p_values.shape)
    pval_corrected.fill(np.nan)
    null_hypothesis.fill(np.nan)
    pval_corrected[mask] = statsmodels.stats.multitest.fdrcorrection(p_values[mask], alpha=0.05, method='indep', is_sorted=False)[1]
    null_hypothesis[mask] = statsmodels.stats.multitest.fdrcorrection(p_values[mask], alpha=0.05, method='indep', is_sorted=False)[0]
    #converts lists with statistical values into panda series and adds these as columns to the merged dataframe
    xsquared_values = pd.Series(xsquared_values)
    p_values = pd.Series(p_values)
    pval_corrected = pd.Series(pval_corrected)
    null_hypothesis = pd.Series(null_hypothesis)
    mdf['FDR_pvalues'] = pval_corrected.values
    mdf['null_hyp_test'] = null_hypothesis.values
    mdf['x_sq_values'] = xsquared_values.values
    mdf['p_values'] = p_values.values
    return mdf
        
# command-line executable script
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--sample1', help = "File of sample 1.")
parser.add_argument('-2', '--sample2', help = "File of sample 2.")
""" parser.add_argument('-o', '--outstem', help = "Stem name for output file. Default stem is 'output' if not given.") """

args = parser.parse_args()

# calls function and returns to variable 'dataframe'
df = stats_dataframe(args.sample1, args.sample2)
df['genes'] = l
print(df)

""" df = pd.DataFrame(np.random.rand(50, 4), columns=['a', 'b', 'c', 'd']) """
import matplotlib.pyplot as plt
""" sampleA = dataframe['sampleA']
sampleB = dataframe['sampleB']
print(l)
ax = plt.scatter(x=l, y=sampleA, color = 'DarkBlue', label = 'Group 1')
plt.scatter(x=l, y=sampleB, color = 'DarkGreen', label = 'Group 2', ax=ax)
plt.show() """

df.plot.scatter(x='sampleA', y='sampleB', color='DarkGreen')
plt.show()
""" To plot multiple column groups in a single axes, repeat plot method specifying target ax. It is recommended to specify color and label keywords to distinguish each groups.
In [62]: ax = df.plot.scatter(x='a', y='b', color='DarkBlue', label='Group 1');
In [63]: df.plot.scatter(x='c', y='d', color='DarkGreen', label='Group 2', ax=ax);
"""