
# NGSaneuploidy

Tool for inferring ploidy levels, testing for aneuploidy and localised CNV.

The following are required:

* python 3, with packages `gzip, sys, numpy, argparse, random, pickle, pathlib, math, scipy.stats, statistics, sklearn`
* R, with packages `pracma, data.table, Rcpp`

Included in this package in addition to runnable scripts are:

*`generics.py` - Contains useful functions used throughout scripts
*`Polynomial_Regressor.pk1` - Fitted variable transformation to process variables ready for aneuploidy test
*`Aneuploidy_Classifier.sav` - Fitted trained logistic regression classifier used in aneuploidy detection


## Generating genotype likelihoods

Overview: calculate genotype likelihoods

`Genotype_Likelihoods.py names.filelist`

### Input

* `Input`: name of a text file containing the suffix of each `.mpileup.gz` file

### Options

* `-o` or `--outFolder`: Output folder. Default: the folder of each input files
* `-i` or `--Inbreeding`: Inbreeding coefficients for each sample accepted as a comma seperated list e.g `0.3,0.2,0.1` alternatively can take in the format `0.2x3,0.4` which is equivilent to `0.2,0.2,0.2,0.4`. All values must be between 0 and 1. Default value is `0xNSAMS`
* `-d` or `--downsampling`: Fraction of the data to be included included in the calculation of genotype likelihoods and aneuploidy inference. That is for a value `v` in [0,1] for each read there is a `vx100%` chance the base is included in the calculations. this can be used to speed up calculations for high coverage samples. Be careful using this argument for low coverage data. Default: `1`
* `-m` or `--min_non_major_freq`: Set the minimum frequency of non major alleles for bases to be included in the calculations. Default: `0.2`
* `-M2` or `--max_minor2_freq`: Set the maximum frequency of third most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `-M3` or `--max_minor3_freq`: Set the maximum frequency of fourth most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `-dp` or `--min_global_depth`: Set the minimum global depth of a base to be included in calculations. All bases with more than this number of reads, after filtering for other conditions mentioned above, across all bases will be included.

### Output

A `.genolikes.gz` file for each prefix in the input file. The columns of the file represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, major allele counts, minor allele counts, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 8 (9 columns)

### Syntax Example

```Shell
python Genotype_Likelihoods.py names.filelist -i 0.1x7,0.15,0.1x2 -d 0.9 -m 0.2 -M2 0.15 -M3 -0.1 -dp 5
```


## Testing for aneuploidy and detecting CNV

Genotype likelihood data can be porcessed and analysed to detect cases of aneuploidy within a set of samples. Running `Aneuploidy_Test.py` will determine a baseline ploidy level for the set of samples and return a list of samples inferred to have aneuploidy with individually inferred ploidy levels for each. 

### Inputs

`Aneuploidy_Test.py` is run with the required arguments of:
* `Input`: name of a text file containing the suffix of each `.mpileup.gz` file used in analysis

Additionally the following options are available:
* `--window [-w]`: The size of the windows to be used in the analysis (number of SNPs per a window)

### Outputs

On screen results will show the probability of aneuploidy before and after each sample is removed until aneuploidy is no longer detected. 

In addition to the on-screen outputs the following output file is produced:


#### `.ploids`
Outputted are resut in the form of a `.ploids` file that contains the following information. First line is the inferred ploidy of each sample. Second line states whether each sample was detected to have aneuploidy or not, 0 refers to no CCNV and 1 means CCNV was inferred. All following lines provide a window by window inferred ploidy for each samples (one line per a sample).

### Syntax Example

```Shell
python Aneuploidy_Test.py names.filelist -w 100
```

## Application to cancer genetics
It has been shown the CNV can promote the production and reoccurances of cancerous tumours. Testing patients all with tumours it is possible to look for genes displaying CNV across all samples for further study. By applying `Aneuploidy_cancer.py` to a genotype likelihoods file it is possible to do this. 

### Inputs

`Aneuploidy_Cancer.py` is run with the required arguments of:
* `Input`: name of a text file containing the suffix of each `.mpileup.gz` file used in analysis

Additionally the following options are available:
* `--window [-w]`: The size of the windows to be used in the analysis (number of SNPs per a window). Choose to approximately match the size of genes of interest

### Outputs

In addition to the on screen outputs displayed from running `Aneuploidy_Test.py`  and production of `.ploids` file (see above) used in visualisations, additional onscreen results show a vector of the starting bases of the windows containing CNV across all samples


### Syntax Example

```Shell
python Aneuploidy_Cancer.py names.filelist -w 100
```

## Visualising results
There are four seperate ways to visulaise the results from the Aneuploidy test. These are:
* `Visualise_multiple.R` where each samples is displayed individually showing the variation in genotype likelihood across the chromosome. All samples are compiled together in the format of a pdf with a page for each plot
* `Visualise_individuals.R` Here tests from multiple chromosomes are compiled to demonstrate any changes in ploidy across the entire genome. Each sample is produces a pdf containing the normalised depths across the chromosomes with the inferred ploidy level for each
* `Visualise_simulations.R` Identical to `Visualise_individuals.R` but optimised for use with simulated data
* `Visualise_cancer.R` Specifically to be used to visualise the results of `Aneuploidy_Cancer.py`. Samples are aligned and localised window ploidies are displayed for a clear visualisation to check for consistent CNV across samples.

### Inputs
Inputs for all methods are similar and are as follows:

#### `Visualise_multiple.R` 
* `input` : name of a text file containing the suffix of each `.mpileup.gz` file used in analysis
* `window` : The size of the windows used in the analysis (see `Anueploidy_Test.py`)

#### `Visualise_individuals.R` 
* `input` : name of a text file containing the suffix of each `.mpileup.gz` file used in analysis
* `NSAMS` : The number of samples included in the analysis
* `window` : The size of the windows used in the analysis (see `Anueploidy_Test.py`)

#### `Visualise_simulations.R` 
* `input` : name of a text file containing the suffix of each `.mpileup.gz` file used in analysis
* `window` : The size of the windows used in the analysis (see `Anueploidy_Test.py`)


#### `Visualise_cancer.R` 
* `input` : name of a text file containing the suffix of each `.mpileup.gz` file used in analysis
* `NSAMS` : The number of samples included in the analysis


### Syntax Example

```Shell
Rscript Visualise_multiple.R names.filelist 100
Rscript Visualise_individuals.R names.filelist 5 100
Rscript Visualise_simulationsR names.filelist 100
Rscript Visualise_cancer.R names.filelist 5
```


## Generating Simulated data
NGS data can be simulated in the form an mpileup file by using `simulMpileu.R`. Not I do not own this script although I have edited it to allow for the inclusion of inbreeding. For full details on how to use the script please visit https://github.com/ImperialCollegeLondon/ngsJulia. 
In addition to the arguments presented I have add the following compulsory argment:
* `--inbreed`: The inbreeding coefficient per sample I.e the proportion of the population that is inbred from which each sample was extracted. Takes input in the same format as `--copy` e.g 0.2x3,0.1 is 0.2,0.2,0.2,0.1". 

### Example 
Simulated data has been provied for practice with the scripts. This can be found saved in the folder `Example` with the corresponding file containing their suffixes callsed `basenames_test` located in the main folder. Each simulation consists of 1 chromosome with 10 samples whose ploidy level is displayed in the fiel name. e.g. T1x4,3,4,1x3,2 consists of 4 haploids, a triploid, a tetraploid, 3 haploids and a diploid. 

### Syntax example of using example data
```Shell
python Genotype_Likelihoods.py basenames_test -o ./Example 
python Aneuploidy_Test.py basenames_test -w 100
Rscript Visulalise_simulations.R basenames_test 100
```

