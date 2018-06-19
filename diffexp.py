#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.ma as ma # for masking functions
import scipy
from scipy.stats import chi2_contingency
import statsmodels.stats.multitest
import sys

def loadf(file) :
    # loads in expression value files, first row are headers
    df = pd.read_csv(file, sep = "\t")
    global l
    l = df['gene'].tolist()
    df.index = l
    del df['gene']
    return df

def merge(df1, df2) :
    # merges two dataframes, if gene row exists in one but not other dataframe, gene will be added with 'NaN' value for sample without gene
    df = df1.join(df2, how='outer')
    return df

# includes both'merge' and 'loadf' function in one
# function which creates pandas data frame containing list of genes and expression values
# and adds columns containing calculated chi x-squared values and p-values

def stats_dataframe(file1, file2) :
    # loads in two expression value files and creates dataframes by calling function 'loadf'
    try :
        f1 = loadf(file1)
        f2 = loadf(file2)
    except :
        print("Error: file loading failed. Ensure two tab-delimited files were entered using appropriate flags. Ensure column headed 'gene' exists containing gene names. See README.md for further information on usage.")
        sys.exit()
    # extracts the two sample labels for index and plot purposes
    global sample_index
    try :
        sample_index = list(f1.columns.values)
        sample_index.append(f2.columns.values[0])
    # checks that both samples were labelled
    except :
        print("Error: Samples may be incorrectly labelled. Ensure both sample files are correctly headed. See README.md for further information on usage.")
        sys.exit()
    # checks that both samples were labelled
    if len(sample_index) != 2 :
        print("Error: Samples are not correctly labelled. Ensure both sample files are correctly headed. See README.md for further information on usage.")
        sys.exit()
    # checks that sample names are not identical
    if sample_index[0] == sample_index[1] :
        print("Error: Sample headers are identical. Ensure sample files have different header sample labels. See README.md for further information on usage.")
        sys.exit()
    # extracts gene labels from each expression file and checks if all gene labels are equal
    f1_genes = sorted(f1.index[0:].format())
    f2_genes = sorted(f2.index[0:].format())
    if f1_genes != f2_genes :
        print("Error: Gene labels are not identical in expression value files loaded.")
        sys.exit()
    # merges the two dataframes into one
    mdf = merge(f1, f2)
    # sums the total of each expession column for statistical calculations
    try:
        A = (sum(mdf[sample_index[0]]))
        B = (sum(mdf[sample_index[1]]))
    except:
        print("Error: Failed to sum expression values. Please ensure values are numerical only.")
        sys.exit()
    # creates empty lists for statistical values to be appended to
    xsquared_values = []
    p_values =[]
    # creates numpy matrix from merged dataframe
    np_df = mdf.as_matrix()
    # iterates over each row of matrix to calculate statistical values for each row
    for i in range(0, len(np_df)) :
        a = np_df[i,0]
        b = np_df[i,1]
        # 'if' 'else' statement which prevents error if two expression values for same gene = zero
        if a == 0 and b == 0 :
            xsquared_values.append(np.nan)
            p_values.append(np.nan)
            print("Warning: Expession values were zero from both samples for gene %s, p-values will be 'NaN' for these and they will not appear on the plot" % (mdf.index[[i]].format()))
        else :
            # calculates x-squared and p-values using chi contigency table and appends values to previously defined lists
            x = np.array([[a,b], [A, B]])
            chi = chi2_contingency(x)
            xsquared_values.append(chi[0])
            p_values.append(chi[1]) 
    # creates numpy array out of p_values from for loop   
    p_values = np.array(p_values)
    # creates Boolean numpy array with p_values with TRUE if not 'NaN'
    mask = np.isfinite(p_values)
    # creates empty numpy arrays with sample???
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
    mdf['x_sq_values'] = xsquared_values.values
    mdf['p_values'] = p_values.values
    mdf['FDR_pvalues'] = pval_corrected.values
    mdf['null_hyp_test'] = null_hypothesis.values
    if args.datafile :
        mdf.to_csv(args.datafile + ".csv", header=True, sep=',', mode='w')
    return mdf

def makeplot(dataframe) :
    colors = np.array(df["null_hyp_test"])
    colors2 = np.empty(colors.shape, dtype=str)
    for c in range(0, len(colors)) :
        if colors[c] == True :
            colors2[c] = 'blue'
        elif colors[c] == False :
            colors2[c] = 'red'
        else :
            colors2[c] = 'black'
    df.plot.scatter(sample_index[0], sample_index[1], c=colors2)  # add colour labels
    if args.outstempdf :
        plt.savefig(args.outstempdf + ".pdf")
    elif args.outstemjpeg :
        plt.savefig(args.outstemjpeg + ".jpeg")
    else:
        print("Error: Please specify output plot format and name, with -opdf or -ojpeg followed by stem name of the plot.")

def mainfunction() :
    global df
    df = stats_dataframe(args.sample1, args.sample2)
    makeplot(df)

# command-line executable script
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--sample1', help = "File for sample 1.")
parser.add_argument('-2', '--sample2', help = "File for sample 2.")
parser.add_argument('-p', '--outstempdf', help = "Stem name for output pdf plot. Must be specified.")
parser.add_argument('-j', '--outstemjpeg', help = "Stem name for output jpeg plot. Must be specified.")
parser.add_argument('-d', '--datafile', help = "Stem name for output csv file. Must be specified.")
args = parser.parse_args()

if args.sample1 and args.sample2 and args.outstempdf :
    # calls function and returns to variable 'df'
    mainfunction()
else :
    print("Error: Necessary arguments have not been provided. Please see README.md for usage information.")