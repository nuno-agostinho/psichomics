# Parse sample information from TCGA sample identifiers

Parse sample information from TCGA sample identifiers

## Usage

``` r
parseTCGAsampleTypes(
  samples,
  filename = system.file("extdata", "TCGAsampleType.RDS", package = "psichomics")
)

parseTCGAsampleInfo(samples, match = NULL)
```

## Arguments

- samples:

  Character: sample identifiers

- filename:

  Character: path to RDS file containing corresponding types

- match:

  Integer: match between samples and subjects (`NULL` by default;
  performs the match)

## Value

Metadata associated with each TCGA sample

## See also

Other functions associated with TCGA data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`getTCGAdataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md),
[`isFirebrowseUp()`](https://nuno-agostinho.github.io/psichomics/reference/isFirebrowseUp.md),
[`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md)

## Examples

``` r
parseTCGAsampleTypes(c("TCGA-01A-Tumour", "TCGA-10B-Normal"))
#>                     01                     10 
#>  "Primary solid Tumor" "Blood Derived Normal" 
samples <- c("TCGA-3C-AAAU-01A-11R-A41B-07", "TCGA-3C-AALI-01A-11R-A41B-07",
             "TCGA-3C-AALJ-01A-31R-A41B-07", "TCGA-3C-AALK-01A-11R-A41B-07",
             "TCGA-4H-AAAK-01A-12R-A41B-07", "TCGA-5L-AAT0-01A-12R-A41B-07")

parseTCGAsampleInfo(samples)
#>                                     Sample types   Patient ID
#> TCGA-3C-AAAU-01A-11R-A41B-07 Primary solid Tumor TCGA-3C-AAAU
#> TCGA-3C-AALI-01A-11R-A41B-07 Primary solid Tumor TCGA-3C-AALI
#> TCGA-3C-AALJ-01A-31R-A41B-07 Primary solid Tumor TCGA-3C-AALJ
#> TCGA-3C-AALK-01A-11R-A41B-07 Primary solid Tumor TCGA-3C-AALK
#> TCGA-4H-AAAK-01A-12R-A41B-07 Primary solid Tumor TCGA-4H-AAAK
#> TCGA-5L-AAT0-01A-12R-A41B-07 Primary solid Tumor TCGA-5L-AAT0
```
