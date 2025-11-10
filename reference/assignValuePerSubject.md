# Assign average sample values to their corresponding subjects

Assign average sample values to their corresponding subjects

## Usage

``` r
assignValuePerSubject(
  data,
  match,
  clinical = NULL,
  patients = NULL,
  samples = NULL
)
```

## Arguments

- data:

  One-row data frame/matrix or vector: values per sample for a single
  gene

- match:

  Matrix: match between samples and subjects

- clinical:

  Data frame or matrix: clinical dataset (only required if the
  `subjects` argument is not handed)

- patients:

  Character: subject identifiers (only required if the `clinical`
  argument is not handed)

- samples:

  Character: samples to use when assigning values per subject (if
  `NULL`, all samples will be used)

## Value

Values per subject

## See also

Other functions to analyse survival:
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
[`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md),
[`survfit.survTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md),
[`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)

## Examples

``` r
# Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
annot <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")

psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 

# Match between subjects and samples
match <- rep(paste("Subject", 1:3), 2)
names(match) <- colnames(psi)

# Assign PSI values to each subject based on the PSI of their samples
assignValuePerSubject(psi[3, ], match)
#> Subject 1 Subject 2 Subject 3 
#> 0.4336962 0.4815878 0.3925439 
```
