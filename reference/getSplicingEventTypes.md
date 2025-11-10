# Get supported splicing event types

Get supported splicing event types

## Usage

``` r
getSplicingEventTypes(psi = NULL, acronymsAsNames = FALSE)
```

## Arguments

- psi:

  Data frame or matrix: alternative splicing quantification data

- acronymsAsNames:

  Boolean: return acronyms as names?

## Value

Named character vector with splicing event types

## See also

Other functions for PSI quantification:
[`filterPSI()`](https://nuno-agostinho.github.io/psichomics/reference/filterPSI.md),
[`listSplicingAnnotations()`](https://nuno-agostinho.github.io/psichomics/reference/listSplicingAnnotations.md),
[`loadAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/loadAnnotation.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md),
[`quantifySplicing()`](https://nuno-agostinho.github.io/psichomics/reference/quantifySplicing.md)

## Examples

``` r
getSplicingEventTypes()
#>                                          Skipped exon 
#>                                                  "SE" 
#>                               Mutually exclusive exon 
#>                                                 "MXE" 
#>                            Alternative 5' splice site 
#>                                                "A5SS" 
#>                            Alternative 3' splice site 
#>                                                "A3SS" 
#>                                Alternative first exon 
#>                                                 "AFE" 
#>                                 Alternative last exon 
#>                                                 "ALE" 
#> Alternative first exon (exon-centred - less reliable) 
#>                                            "AFE_exon" 
#>  Alternative last exon (exon-centred - less reliable) 
#>                                            "ALE_exon" 
```
