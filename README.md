# diffexp
A programme which determines differential gene expression between a list of genes from two samples, outputting a scatter plot and statistical values. Takes two .txt file containing same genes with read counts, and outputs a pdf scatter plot of read counts against each other from both samples, with colour of the plotted datapoint indicating whether the expression difference is significant. Determines significance using the null hypothesis test and false discovery rate-corrected p-values. A table containing statistical values is outputted as a csv file.

## Dependencies

1. Python3

2. The following Python packages:
  * matplotlib
  * pandas
  * numpy
  * scipy
  * statsmodels


## File inputs

* Two tab-delimited text files containing expression values (required): A header row must contain 'gene' for the gene identifier column and the sample name for the second column containing expression values. See example files 'sampleA' and 'sampleB'. Expression values must be numerical. Both files must contain the same genes labels, and so the same number of rows by default.


## File outputs

* PDF Plot : A scatter plot. Each datapoint is expression value of one sample against other sample. Red datapoints represent significant expression differences (null-hypothesis = FALSE) according to an FDR (false-discovery-rate) p-value =< 0.05. Blue datapoints represent no significant difference (null-hypothesis = TRUE).


* CSV File (optional): A script-assembled table containing input genes and expression values, as well as, for each gene, X-squared values, p-values, FDR-corrected p-values and the null-hypothesis boolean value (0=TRUE or 1=FALSE).


## Usage
Two sample files must be given using '-1' and '-2' flag. '-p' flag is required, used to specify the plot stem name. '-d' flag is not required, used to specify the csv stem name.

    python ./diffexp.py -1 *sample1.txt* -2 *sample2.txt* -p *plotname* -d *csvname*


## Statistical Analysis
p-values are calculated using a 2x2 chi-square test using the 'scipy' package 'stats.chi2_contingency' function. To correct for multiple testing, false discover rate (FDR) corrected p-values are calculated. FDR-corrected p-values are calculated from the p-values (disregarding NaN values for genes with reads of '0' from both samples) using the 'statsmodels' package 'stats.multitest' function. The null-hypothesis test uses an alpha=0.05. Boolean values are returned. '0' = null hypothesis is TRUE. '1' = null hypothesis is FALSE.
   
## Example
Example test files are available in the directory 'ExampleFiles': 'sampleA.txt' and 'sampleB.txt' and 'sampleBprime.txt'.

To run script to generate pdf plot and csv file (must be in working directory containing diffexp.py):

    python ./diffexp.py -1 sampleA.txt -2 sample2.txt -p plot -d datatable

Running this script creates 'plot.pdf' and 'datatable.csv', adding them to the working directory. A warning will appear: 'Warning: Expression values were zero from both samples for gene ['EfM3.000130'], p-values will be 'NaN' for these and they will not appear on the plot". This indicates that the gene with expression values of '0' for each sample are not included in the analysis.

The example file 'sampleBprime' contains the following:

* A gene which is not found in sampleA.txt ('EfM3.blimey')
* Two read values which are not in numerical form ('1e8' and 'bob')

Running the above code with sampleBprime.txt in place of sampleB.txt will call an error message 'Error: Gene labels are not identical in expression value files loaded'. The order of the genes in sampleBprime.txt is not identical to sampleA.txt (e.g. EfM3.000010 is at the bottom and not the top of the list), but this will not effect the program as the script does not require identical genes to be in the same order.