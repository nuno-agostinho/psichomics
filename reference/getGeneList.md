# Get curated, literature-based gene lists

Available gene lists:

- **Sebestyen et al., 2016**: 1350 genes encoding RNA-binding proteins,
  167 of which are splicing factors

## Usage

``` r
getGeneList(genes = NULL)
```

## Arguments

- genes:

  Vector of characters: intersect lists with given genes (lists with no
  matching genes will not be returned)

## Value

List of genes

## See also

Other functions for data grouping:
[`createGroupByAttribute()`](https://nuno-agostinho.github.io/psichomics/reference/createGroupByAttribute.md),
[`getSampleFromSubject()`](https://nuno-agostinho.github.io/psichomics/reference/getSampleFromSubject.md),
[`getSubjectFromSample()`](https://nuno-agostinho.github.io/psichomics/reference/getSubjectFromSample.md),
[`groupPerElem()`](https://nuno-agostinho.github.io/psichomics/reference/groupPerElem.md),
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md),
[`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)

## Examples

``` r
getGeneList()
#> Sebestyen et al. 2016
#>   -> Human RNA-binding protein splicing factors [167 genes]: A1CF, ANKHD1, CELF1, CELF2, ...
#>   -> Human RNA-binding proteins [1350 genes]: A1CF, AATF, ABCF1, ABT1, ...
#> 
#> Source: 
#> Endre Sebestyén et al. (2016). Large-scale analysis of genome and transcriptome alterations in multiple tumors unveils novel cancer-relevant splicing networks. Genome Research, 26(6), 732-744
#> ================================================================================
```
