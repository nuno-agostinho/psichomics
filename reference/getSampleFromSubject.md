# Get samples matching the given subjects

Get samples matching the given subjects

## Usage

``` r
getSampleFromSubject(
  patients,
  samples,
  clinical = NULL,
  rm.NA = TRUE,
  match = NULL,
  showMatch = FALSE
)
```

## Arguments

- patients:

  Character or list of characters: subject identifiers

- samples:

  Character: sample identifiers

- clinical:

  Data frame or matrix: clinical dataset

- rm.NA:

  Boolean: remove missing values?

- match:

  Integer: vector of subject index with the sample identifiers as name
  to save time (optional)

- showMatch:

  Boolean: show matching subject index?

## Value

Names of the matching samples (if `showMatch = TRUE`, a character with
the subjects as values and their respective samples as names is
returned)

## See also

Other functions for data grouping:
[`createGroupByAttribute()`](https://nuno-agostinho.github.io/psichomics/reference/createGroupByAttribute.md),
[`getGeneList()`](https://nuno-agostinho.github.io/psichomics/reference/getGeneList.md),
[`getSubjectFromSample()`](https://nuno-agostinho.github.io/psichomics/reference/getSubjectFromSample.md),
[`groupPerElem()`](https://nuno-agostinho.github.io/psichomics/reference/groupPerElem.md),
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md),
[`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)

## Examples

``` r
subjects <- c("GTEX-ABC", "GTEX-DEF", "GTEX-GHI", "GTEX-JKL", "GTEX-MNO")
samples <- paste0(subjects, "-sample")
clinical <- data.frame(samples=samples)
rownames(clinical) <- subjects
getSampleFromSubject(subjects[c(1, 4)], samples, clinical)
#> [1] "GTEX-ABC-sample" "GTEX-JKL-sample"
```
