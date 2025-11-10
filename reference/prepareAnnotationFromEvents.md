# Prepare annotation from alternative splicing events

In case more than one data frame with alternative splicing events is
given, the events are cross-referenced according to the chromosome,
strand and relevant coordinates per event type (see details).

## Usage

``` r
prepareAnnotationFromEvents(...)
```

## Arguments

- ...:

  Data frame(s) of alternative splicing events to include in the
  annotation

## Value

List of data frames with the annotation from different data frames
joined by event type

## Details

Events from two or more data frames are cross-referenced based on each
event's chromosome, strand and specific coordinates relevant for each
event type:

- Skipped exon: constitutive exon 1 end, alternative exon (start and
  end) and constitutive exon 2 start

- Mutually exclusive exon: constitutive exon 1 end, alternative exon 1
  and 2 (start and end) and constitutive exon 2 start

- Alternative 5' splice site: constitutive exon 1 end, alternative exon
  1 end and constitutive exon 2 start

- Alternative first exon: same as alternative 5' splice site

- Alternative 3' splice site: constitutive exon 1 end, alternative exon
  1 start and constitutive exon 2 start

- Alternative last exon: same as alternative 3' splice site

## Note

When cross-referencing events, gene information is discarded.

## See also

Other functions to prepare alternative splicing annotations:
[`parseSuppaAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md)

## Examples

``` r
# Load sample files (SUPPA annotation)
folder <- "extdata/eventsAnnotSample/suppa_output/suppaEvents"
suppaOutput <- system.file(folder, package="psichomics")

# Parse and prepare SUPPA annotation
suppa <- parseSuppaAnnotation(suppaOutput)
#> 
#> Retrieving SUPPA annotation...
#> 
#> Parsing SUPPA annotation...
annot <- prepareAnnotationFromEvents(suppa)
#> 
#> Sorting coordinates...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==========                                                            |  14%  |                                                                              |====================                                                  |  29%  |                                                                              |==============================                                        |  43%  |                                                                              |========================================                              |  57%  |                                                                              |==================================================                    |  71%  |                                                                              |============================================================          |  86%  |                                                                              |======================================================================| 100%
#> Joining events per event type...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==========                                                            |  14%  |                                                                              |====================                                                  |  29%  |                                                                              |==============================                                        |  43%  |                                                                              |========================================                              |  57%  |                                                                              |==================================================                    |  71%  |                                                                              |============================================================          |  86%  |                                                                              |======================================================================| 100%
#> Cleaning the annotation...

# Load sample files (rMATS annotation)
folder <- "extdata/eventsAnnotSample/mats_output/ASEvents/"
matsOutput <- system.file(folder, package="psichomics")

# Parse rMATS annotation and prepare combined annotation from rMATS and SUPPA
mats <- parseMatsAnnotation(matsOutput)
#> 
#> Retrieving rMATS annotation...
#> 
#> Parsing rMATS annotation...
annot <- prepareAnnotationFromEvents(suppa, mats)
#> 
#> Sorting coordinates...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |==========                                                            |  14%  |                                                                              |===============                                                       |  21%  |                                                                              |====================                                                  |  29%  |                                                                              |=========================                                             |  36%  |                                                                              |==============================                                        |  43%  |                                                                              |===================================                                   |  50%  |                                                                              |========================================                              |  57%  |                                                                              |=============================================                         |  64%  |                                                                              |==================================================                    |  71%  |                                                                              |=======================================================               |  79%  |                                                                              |============================================================          |  86%  |                                                                              |======================================================================| 100%
#> Joining events per event type...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==========                                                            |  14%  |                                                                              |====================                                                  |  29%  |                                                                              |==============================                                        |  43%  |                                                                              |========================================                              |  57%  |                                                                              |==================================================                    |  71%  |                                                                              |============================================================          |  86%  |                                                                              |======================================================================| 100%
#> Cleaning the annotation...
```
