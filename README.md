# psichomics [![Build Status](https://travis-ci.com/nuno-agostinho/psichomics.svg?token=WnQvrH4wCa4UkqjJquSq&branch=master)](https://travis-ci.com/nuno-agostinho/psichomics) [![codecov](https://codecov.io/gh/nuno-agostinho/psichomics/branch/master/graph/badge.svg?token=huZOun5jSD&branch=master)](https://codecov.io/gh/nuno-agostinho/psichomics)
R package to analyse and visualise alternative splicing data

## Start running
If you want to install from [Bioconductor](https://www.bioconductor.org), you're out of luck. This package is not 
available in Bioconductor (yet?).

To start using this program, simply type the following in RStudio or in a R console:
```r
library(shiny)
runGitHub("psichomics", "nuno-agostinho")
```

You can also download it and install it as a package. Or clone this repository, set the downloaded folder as the working
directory and type `psichomics()`.

## Data input
### Downloading TCGA data
You can download data from The Cancer Genome Atlas (TCGA) by using the Firehose API incorporated in the package. Simply
choose the cohort of interest, date of the sample, type of interest and so on. Wait for the data download and process to
finish.

### Loading user files
To load your own files, simply choose the folder where you have your data. Psichomics will try to process all the data
contained in the indicated folder and sub-folders.

## Data groups
To create groups for a dataset, first select the category and the tab of the dataset of interest. Next, select how to
create the group:

- **By column:** automatically create groups by selecting a specific column of the dataset; for example, with the `Clinical data`
dataset selected, start typing `pathologic stage`, select the appropriate field from the suggestions, click `Create group` and 
confirm that there is now  one group for each pathologic stage.
- **By row:** input specific rows to create a group
- **By subset expression:** type a subset expression
- **By GREP expression:** apply a GREP expression over a specific column of the dataset

You can also select groups using the checkbox to merge, intersect or remove them.

When changing the dataset, the interface will update accordingly. Worry not, everything is saved.

## Exon/intron inclusion levels
This is also known as percentage spliced-in, percent of splicing index, PSI or Ψ. It is the ratio of read counts for
the isoforms including a specific exon over the read counts of the isoforms that both include and exclude that exon.

If you provide junction read counts, this can be quickly calculated. You can also load a dataset with this information by
using more accurate software such as [MISO](http://genes.mit.edu/burgelab/miso), [VAST-TOOLS](https://github.com/vastgroup/vast-tools),
[rMATS](http://rnaseq-mats.sourceforge.net) or [SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa).

## Author
Nuno Agostinho, 2015 @ <a href="http://imm.medicina.ulisboa.pt/group/compbio/" target="_blank">Nuno Morais Lab, IMM</a>

Special thanks to my lab colleagues for their work-related support and supporting chatter.
