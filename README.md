# psichomics [![Build Status](https://travis-ci.com/nuno-agostinho/psichomics.svg?token=WnQvrH4wCa4UkqjJquSq&branch=master)](https://travis-ci.com/nuno-agostinho/psichomics) [![codecov](https://codecov.io/gh/nuno-agostinho/psichomics/branch/master/graph/badge.svg?token=huZOun5jSD)](https://codecov.io/gh/nuno-agostinho/psichomics)
R package to analyse and visualise alternative splicing data

## Start running
If you want to install from [Bioconductor](https://www.bioconductor.org), you're
out of luck. This package is not available in Bioconductor (yet?).

To start using this program, follow these steps:

1. [Install R](https://www.r-project.org/)
2. Open [RStudio](https://www.rstudio.com/products/rstudio) (or open a console, 
type `R` and press enter)
3. Type the following:
```r
install.packages("devtools")
devtools::install_github("nuno-agostinho/psichomics", auth_token = "b752621b2059cd64ac21d5f7e8418821feb81b81")
library(psichomics)
psichomics()
```

You can also download it and install it as a package. Or clone this repository, 
set the downloaded folder as the working directory and type
`devtools::load.all(); psichomics()`.

## Data input
### Downloading TCGA data
You can download data from
[The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov) by using the
Firehose API incorporated in the package. Simply choose the cohort of interest, 
date of the sample, type of interest and so on. Wait for the downloads to finish
and then click again in "Get Data" to process and load the data.

### Loading user files
To load your own files, simply choose the folder where you have your data. 
Psichomics will try to process all the data contained in the indicated folder 
and sub-folders to search for files that can be loaded.

## Exon/intron inclusion levels
This is also known as percentage spliced-in, percent of splicing index, PSI or
Ψ. It is the proportion of reads supporting the inclusion of an exon over the
reads supporting both the inclusion and exclusion of that exon.

If you provide junction read counts, this can be quickly calculated. You can 
also load a dataset with this information by using more accurate software such 
as [MISO](http://genes.mit.edu/burgelab/miso),
[VAST-TOOLS](https://github.com/vastgroup/vast-tools), 
[rMATS](http://rnaseq-mats.sourceforge.net) or 
[SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa).

## Analyses
The program performs survival and principal component analyses, as well as
differential splicing analysis using non-parametric statistical tests.

## Data groups
- **By column:** automatically create groups by selecting a specific column of 
the dataset; for instance, to create a group for each tumour stage, start typing
`tumor_stage`, select the appropriate field from the suggestions, click 
`Create group` and confirm that there is now one group for each stage.
- **By row:** input specific rows to create a group
- **By subset expression:** type a subset expression
- **By GREP expression:** apply a GREP expression over a specific column of the 
dataset

You can also select groups using the checkbox in order to merge, intersect or 
remove the groups.

## Author
Nuno Agostinho, 2015-2016 @ <a href="http://imm.medicina.ulisboa.pt/group/compbio/" target="_blank">Nuno Morais Lab, IMM</a>

Special thanks to my lab colleagues for their work-related support and 
supporting chatter.
