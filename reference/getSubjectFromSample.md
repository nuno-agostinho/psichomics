# Get subjects from given samples

Get subjects from given samples

## Usage

``` r
getSubjectFromSample(sampleId, patientId = NULL, na = FALSE, sampleInfo = NULL)
```

## Arguments

- sampleId:

  Character: sample identifiers

- patientId:

  Character: subject identifiers to filter by (optional; if a matrix or
  data frame is given, its rownames will be used to infer the subject
  identifiers)

- na:

  Boolean: return `NA` for samples with no matching subjects

- sampleInfo:

  Data frame or matrix: sample information containing the sample
  identifiers as rownames and a column named "Subject ID" with the
  respective subject identifiers

## Value

Character: subject identifiers corresponding to the given samples

## See also

Other functions for data grouping:
[`createGroupByAttribute()`](https://nuno-agostinho.github.io/psichomics/reference/createGroupByAttribute.md),
[`getGeneList()`](https://nuno-agostinho.github.io/psichomics/reference/getGeneList.md),
[`getSampleFromSubject()`](https://nuno-agostinho.github.io/psichomics/reference/getSampleFromSubject.md),
[`groupPerElem()`](https://nuno-agostinho.github.io/psichomics/reference/groupPerElem.md),
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md),
[`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)

## Examples

``` r
samples <- paste0("GTEX-", c("ABC", "DEF", "GHI", "JKL", "MNO"), "-sample")
getSubjectFromSample(samples)
#> GTEX-ABC-sample GTEX-DEF-sample GTEX-GHI-sample GTEX-JKL-sample GTEX-MNO-sample 
#>      "GTEX-ABC"      "GTEX-DEF"      "GTEX-GHI"      "GTEX-JKL"      "GTEX-MNO" 

# Filter returned samples based on available subjects
subjects <- paste0("GTEX-", c("DEF", "MNO"))
getSubjectFromSample(samples, subjects)
#> GTEX-DEF-sample GTEX-MNO-sample 
#>      "GTEX-DEF"      "GTEX-MNO" 
```
