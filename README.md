# PSIchomics [![Build Status](https://travis-ci.com/nuno-agostinho/psichomics.svg?token=WnQvrH4wCa4UkqjJquSq&branch=master)](https://travis-ci.com/nuno-agostinho/psichomics) [![codecov](https://codecov.io/gh/nuno-agostinho/psichomics/branch/master/graph/badge.svg?token=huZOun5jSD)](https://codecov.io/gh/nuno-agostinho/psichomics)
R package to quantify, analyse and visualise alternative splicing data from 
[The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov/docs/publications/tcga).

## Available tutorials

The following tutorial is available:

* [Visual interface](http://rpubs.com/nuno-agostinho/psichomics-tutorial-visual)

Other tutorials coming soon:

* Command-line interface
* Developers and other contributors

## Install and start running
This package is not available in [Bioconductor](https://www.bioconductor.org)
yet. To install and start using this program, follow these steps:

1. [Install R](https://www.r-project.org/)
2. Open [RStudio](https://www.rstudio.com/products/rstudio) (or open a console, 
type `R` and press enter)
3. Type the following to install, load and start the visual interface:
```r
install.packages("devtools")
devtools::install_github("nuno-agostinho/psichomics")
library(psichomics)
psichomics()
```

## Data input
### Downloading TCGA data
You can download data from
[The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov) using this
package. Simply choose the cohort of interest, date of the sample, type of 
interest and so on. Wait for the downloads to finish and then click again in 
**Load Data** to process and load the data.

### Loading user files
To load your own files, simply choose the folder where the data is located. 
PSIchomics will try to process all the data contained in the given folder and
sub-folders to search for files that can be loaded.

## Exon/intron inclusion levels
To quantify alternative splicing based on the porportion of isoforms that 
include an exon, the Percent Spliced-In (PSI or Ψ) metric is used.

An estimate of this value is obtained based on the the proportion of reads 
supporting the inclusion of an exon over the reads supporting both the inclusion
and exclusion of that exon. To measure this estimate, both alternative splicing 
annotation and junction quantification are required. While alternative splicing 
Human (hg19 assembly) annotation is already provided, junction quantification 
may be retrieved from TCGA.

## Analyses
The program performs survival and principal component analyses, as well as
differential splicing analysis using non-parametric statistical tests.

### Differential splicing analysis
Analyse alternative splicing quantification based on variance and median 
statistical tests. The groups available for differential analysis comprise 
sample types (e.g. normal versus tumour) and clinical attributes of patients 
(e.g. tumour stage).

### Gene, transcript and protein information
For a given splicing event, examine its gene's annotation and corresponding 
transcripts and proteins. Related research articles are also available.

### Principal component analysis (PCA)
Explore alternative splicing quantification groups using associated clinical 
attributes.

### Survival analysis
Analyse survival based on clinical attributes (e.g. tumour stage, gender and
race). Additionally, study the impact of the quantification of a single 
alternative splicing event on patient survivability.

## Data groups

- **By column:** automatically create groups by selecting a specific column of 
the dataset; for instance, to create a group for each tumour stage, start typing
`tumor_stage`, select the appropriate field from the suggestions, click 
`Create group` and confirm that there is now one group for each stage.
- **By row:** input specific rows to create a group
- **By subset expression:** type a subset expression
- **By GREP expression:** apply a GREP expression over a specific column of the 
dataset

You can also select groups by clicking on them in order to merge, intersect or 
remove the groups.

## Feedback
All feedback on the program, documentation and associated material is welcome. 
Please, send any suggestions and comments to the following contact:

Nuno Saraiva-Agostinho (nunodanielagostinho@gmail.com)
<a href="http://imm.medicina.ulisboa.pt/group/compbio/" target="_blank">Nuno Morais Lab, IMM</a>

Special thanks to my lab colleagues for their work-related support and 
supporting chatter.