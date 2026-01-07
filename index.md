# psichomics

> **Original article:**
>
> Nuno Saraiva-Agostinho and Nuno L. Barbosa-Morais (2019). [psichomics:
> graphical application for alternative splicing quantification and
> analysis](https://doi.org/10.1093/nar/gky888). *Nucleic Acids
> Research*. 47(2), e7.

Interactive R package with an intuitive Shiny-based graphical interface
for alternative splicing quantification and integrative analyses of
alternative splicing and gene expression based on [The Cancer Genome
Atlas (TCGA)](https://tcga-data.nci.nih.gov), the [Genotype-Tissue
Expression (GTEx) project](http://www.gtexportal.org), [Sequence Read
Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) and user-provided data.

*psichomics* interactively performs: - Dimensionality reduction -
Median- and variance-based differential splicing and gene expression
analyses - Survival analysis - Correlation analysis - Grouping by
clinical and molecular features (such as tumour stage or survival) -
Genomic mapping and functional annotation of alternative splicing events
and genes

![Differential splicing analysis in
psichomics](reference/figures/screenshot.png)

Differential splicing analysis in *psichomics*

## Table of Contents

- [Install and start running](#install-and-start-running)
  - [Bioconductor](#bioconductor)
  - [GitHub](#github)
  - [Docker](#docker)
- [Tutorials](#tutorials)
- [Workflow](#workflow)
  - [Data input](#data-input)
    - [Alternative splicing
      quantification](#alternative-splicing-quantification)
    - [Gene expression processing](#gene-expression-processing)
  - [Data grouping](#data-grouping)
  - [Data analyses](#data-analyses)
- [Feedback and support](#feedback-and-support)
- [References](#references)

## Install and start running

### Bioconductor

To install the package from
[Bioconductor](https://www.bioconductor.org), type the following in
[RStudio](https://www.rstudio.com/products/rstudio) or in an R console:

``` r
install.packages("BiocManager")
BiocManager::install("psichomics")
library("psichomics")
```

Start the visual interface of psichomics with
[`psichomics()`](https://nuno-agostinho.github.io/psichomics/reference/psichomics.md)

### GitHub

Install from GitHub (specify a branch or tag via the `ref` argument):

``` r
install.packages("remotes")
remotes::install_github("nuno-agostinho/psichomics", ref="master")
library("psichomics")
```

Start the visual interface of psichomics with
[`psichomics()`](https://nuno-agostinho.github.io/psichomics/reference/psichomics.md)

### Docker

The Docker images are based on [Bioconductor
Docker](https://www.bioconductor.org/help/docker/) and contain
psichomics and its dependencies.

1.  Pull the latest Docker image:

&nbsp;

    docker pull ghcr.io/nuno-agostinho/psichomics:latest

2.  Start RStudio Web from the Docker image (mount your own Downloads
    folder):

&nbsp;

    docker run -e PASSWORD=bioc -p 8787:8787 -v ~/Downloads:/home/rstudio/Downloads ghcr.io/nuno-agostinho/psichomics:latest

3.  Go to RStudio Web via the web browser at <https://localhost:8787>
4.  Log in RStudio with user `rstudio` and password `bioc`
5.  Load psichomics using
    [`library(psichomics)`](https://nuno-agostinho.github.io/psichomics/)
6.  Start the visual interface of psichomics with
    [`psichomics()`](https://nuno-agostinho.github.io/psichomics/reference/psichomics.md)

## Tutorials

The following case studies and tutorials are available and were based on
our [original article](https://doi.org/10.1093/nar/gky888):

- [Visual
  interface](https://nuno-agostinho.github.io/psichomics/articles/GUI_tutorial.html)
- [Command-line
  interface](https://nuno-agostinho.github.io/psichomics/articles/CLI_tutorial.html)
- [Loading user-provided
  data](https://nuno-agostinho.github.io/psichomics/articles/custom_data.html)
- [Preparing alternative splicing
  annotations](https://nuno-agostinho.github.io/psichomics/articles/AS_events_preparation.html)

Another tutorial was published as part of the Methods in Molecular
Biology book series (the code for performing the analysis can be found
[here](https://github.com/nuno-agostinho/stem-cell-analysis-in-psichomics)):

> Nuno Saraiva-Agostinho and Nuno L. Barbosa-Morais (2020).
> **[Interactive Alternative Splicing Analysis of Human Stem Cells Using
> psichomics](https://doi.org/10.1007/978-1-0716-0301-7_10)**. In:
> Kidder B. (eds) Stem Cell Transcriptional Networks. *Methods in
> Molecular Biology*, vol 2117. Humana, New York, NY

## Workflow

### Data input

Automatic retrieval and loading of pre-processed data from the following
sources:

- [TCGA](https://tcga-data.nci.nih.gov) data of given tumours, including
  subject- and sample-associated information, junction quantification
  and gene expression data
- [GTEx](http://www.gtexportal.org) data of given tissues, including
  subject- and sample-associated information, junction quantification
  and gene expression data
- [SRA](https://www.ncbi.nlm.nih.gov/sra) data from select SRA projects
  via the [recount](https://jhubiostatistics.shinyapps.io/recount/)
  package

Other SRA, [VAST-TOOLS](https://github.com/vastgroup/vast-tools) and
user-provided data can also be manually loaded. Please read [Loading
user-provided
data](https://nuno-agostinho.github.io/psichomics/articles/custom_data.html)
for more information.

#### Alternative splicing quantification

The quantification of each alternative splicing event is based on the
proportion of junction reads that support the inclusion isoform, known
as percent spliced-in or PSI [(Wang *et al.*,
2008)](http://www.nature.com/nature/journal/v456/n7221/full/nature07509.md).

An estimate of this value is obtained based on the the proportion of
reads supporting the inclusion of an exon over the reads supporting both
the inclusion and exclusion of that exon. To measure this estimate, we
require:

1.  **Alternative splicing annotation**: human annotation is provided
    and custom annotations can be prepared for use in psichomics.
2.  Quantification of RNA-Seq reads aligning to exon-exon splice
    junctions (**exon-exon junction quantification**), either
    user-provided or retrieved from
    [TCGA](https://tcga-data.nci.nih.gov),
    [GTEx](http://www.gtexportal.org) and
    [SRA](https://www.ncbi.nlm.nih.gov/sra).

#### Gene expression processing

Gene expression can be normalised, filtered and log2-transformed in-app
or provided by the user.

### Data grouping

Molecular and clinical sample-associated attributes allow to establish
groups that can be explored in data analyses.

For instance, [TCGA](https://tcga-data.nci.nih.gov) data can be analysed
based on smoking history, gender and race, among other attributes.
Groups can also be manipulated (e.g. merged, intersected, etc.),
allowing for complex attribute combinations. Groups can also be saved
and loaded between different sessions.

### Data Analyses

- **Dimensionality reduction** via principal and independent component
  analysis (PCA and ICA) on alternative splicing quantification and gene
  expression.

- **Differential splicing and gene expression analysis** based on
  variance and median parametric and non-parametric statistical tests.

- **Correlation between gene expression and splicing quantification**,
  useful to correlate the expression of a given event with the
  expression of RNA-binding proteins, for instance.

- **Survival analysis** via Kaplan-Meier curves and Cox models based on
  sample-associated features. Additionally, we can study the impact of a
  splicing event (based on its quantification) or a gene (based on its
  expression) on patient survivability.

- **Gene, transcript and protein annotation**, including relevant
  research articles.

## Feedback and support

Please send any feedback and questions on psichomics to:

> Nuno Saraiva-Agostinho (<nunodanielagostinho@gmail.com>)
>
> [Disease Transcriptomics Lab, Instituto de Medicina Molecular
> (Portugal)](http://imm.medicina.ulisboa.pt/group/distrans/)

## References

Wang, E. T., R. Sandberg, S. Luo, I. Khrebtukova, L. Zhang, C. Mayr, S.
F. Kingsmore, G. P. Schroth, and C. B. Burge. 2008. [*Alternative
isoform regulation in human tissue
transcriptomes.*](http://www.nature.com/nature/journal/v456/n7221/full/nature07509.md)
Nature 456 (7221): 470–76.
