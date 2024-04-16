# Comparing chromatin contact maps at scale: methods and insights

Preprint link: https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1

In this repo, you can find:
- [Code](https://github.com/pollardlab/contact_map_scoring/blob/main/code/scoring.py) and [Tutorials](https://github.com/pollardlab/contact_map_scoring/tree/main/notebooks) for running scoring methods on both experimental and predicted contact frequency maps.
- [Dataset](https://github.com/pollardlab/contact_map_scoring/tree/main/data) of:
	* Scores for _in silico_ insertions and deletions throughout the genome
	* Scores compartin windows around differentially expressed genes between ESC and HFF in chromosomes 21 and 22

![](Fig1.png)


## Running scoring functions


Download this repo:
```
git clone https://github.com/pollardlab/contact_map_scoring.git
```

We provide scripts to run all 25 scoring functions on experimental maps. Since they vary in input type and coding language, they are split into multiple scripts as follows:


### Methods that use contact maps

Methods that take in matrices that correspond to contact frequency maps of a certain region. This includes the following 13 methods:
- Correlation
- MSE
- SSIM
- Contact directionality (corr)
- Distance enrichment (corr)
- Eigenvector (corr)
- Insulation (corr)
- Insulation (mse)
- Loops
- SCC
- TADs
- Triangle (corr)
- Triangle (mse)

#### Script
[code/run_methods_that_use_contact_maps.py](https://github.com/pollardlab/contact_map_scoring/blob/main/code/run_methods_that_use_contact_maps.py)


#### Installation requirements
- scipy: for Correlation, Contact directionality (corr), Distance enrichment (corr), Eigenvector (corr), Insulation (corr), SCC, and Triangle (corr)
- skimage: for Contact directionality (corr), SSIM, Triangle (corr), and Triangle (mse)
- sklearn: for Eigenvector only
- cooltools: for TADs only
- hicrep: for SCC only


#### Directions
- Generate input files outlined below
- Change variables in the script following instructions there
- Run script in the terminal
```
python run_methods_that_use_contact_maps.py
````

To run these methods on predicted maps instead of experimental ones, follow instructions in [Running_Scoring_Functions_on_Predicted_Maps.ipynb](https://github.com/pollardlab/contact_map_scoring/tree/main/notebooks/Running_Scoring_Functions_on_Predicted_Maps.ipynb)


#### Input
1. Windows to score. This should be a tab-delimited text file with columns: chrom, start, end.
2. Two cool files that will be compared at provided windows.


#### Output 
Tab-delimited text file with the same rows as the input windows file and added columns with scores for each method.



### TADcompare and HiCcompare

#### Script
[code/run_tad_hic_compare.R](https://github.com/pollardlab/contact_map_scoring/blob/main/code/run_tad_hic_compare.R)


#### Installation requirements
- [TADcompare](https://bioconductor.org/packages/release/bioc/vignettes/TADCompare/inst/doc/TADCompare.html#installation)
- [HiCcompare](https://github.com/dozmorovlab/HiCcompare?tab=readme-ov-file#installation)


#### Directions
- Generate input files following instructions below
- Change variables in the script following instructions there
- Run script in the terminal
```
Rscript run_tad_hic_compare.R
````

#### Input
Two chromosome-specific text files generated from cool files. Input file format: {input_file_prefix}_{experiment}_{chromosome}_{resolution}.txt

Example of generating chromosome-specific text file for:
- input_file_prefix: example_cool
- experiment: MicroC
- chromosome: chr21
- resolution: 2048

1. Create chromosome-wide text file from cool file:
```
cooler dump --join example_cool_MicroC.mcool::/resolutions/2048 > example_cool_MicroC_2048.txt
```

2. Create chromosome-specific text file from chromosome-wide text file:
```
awk -F '\t' '$1 == "chr21"{ print }' example_cool_MicroC_2048.txt > example_cool_MicroC_chr21_2048.txt
```


#### Output 
Tab-delimited text files with results for each method, experiment, chromosome, and resolution saved in a new directory created separately for HiC compare and TADcompare



## Contact

Please contact raise an issue on this repo or email us `{laura.gunsalus, evonne.mcarthur}@gmail.com`) if you have any questions or concerns.
