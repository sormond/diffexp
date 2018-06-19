# diffexp
A programme to determine differential expression of genes between two samples

## Dependencies
- Python 3
- The following Python packages:
	matplotlib.pyplot
	pandas
	numpy
	scipy
	statsmodels


## File inputs

* Two tab-delimited text files containing expression values (required): A header row must contain 'gene' for the gene identifier column and the sample name for the second column containing expression values. See example files 'sampleA' and 'sampleB'. Expression values must be numerical. Both files must contain the same genes labels, and so the same number of rows by default.


## File outputs

* PDF Plot : A scatter plot. Each datapoint is expression value of one sample against other sample. Blue datapoints represent significant expression differences (null-hypothesis = FALSE) according to an FDR p-value < 0.05. Red datapoints represent no significant difference (null-hypothesis = TRUE).


* CSV File (optional): A script-assembled table containing input genes and expression values, as well as, for each gene, X-squared values, p-values, FDR-corrected p-values and the null-hypothesis boolean value (0=TRUE or 1=FALSE).



## Usage
Two sample files must be given using '-1' and '-2' flag. '-p' flag is required, used to specify the plot stem name. '-d' flag is not required, used to specify the csv stem name.

   ./diffexp.py -1 *sample1.txt* -2 *sample2.txt* -p *plotname* -d *csvname*

   
## Example
Example test files are avialable in the directory 'ExampleFiles': 'sampleA.txt' and 'sampleB.txt'.

To run script to generate pdf plot and csv file (must be in working directory containing diffexp.py):

   ./diffexp.py -1 sampleA.txt -2 sample2.txt -p plot -d datatable

Running this script creates 'plot.pdf' and 'datatable.csv' and adds to the working directory. A warning will appear: 'Warning: Expession values were zero from both samples for gene ['EfM3.000130'], p-values will be 'NaN' for these and they will not appear on the plot".