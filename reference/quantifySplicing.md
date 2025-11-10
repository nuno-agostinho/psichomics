# Quantify alternative splicing events

Quantify alternative splicing events

## Usage

``` r
quantifySplicing(
  annotation,
  junctionQuant,
  eventType = c("SE", "MXE", "ALE", "AFE", "A3SS", "A5SS"),
  minReads = 10,
  genes = NULL
)
```

## Arguments

- annotation:

  List of data frames: annotation for each alternative splicing event
  type

- junctionQuant:

  Data frame: junction quantification

- eventType:

  Character: splicing event types to quantify

- minReads:

  Integer: values whose number of total supporting read counts is below
  `minReads` are returned as `NA`

- genes:

  Character: gene symbols for which to quantify splicing events (if
  `NULL`, events from all genes are quantified)

## Value

Data frame with the quantification of the alternative splicing events

## See also

Other functions for PSI quantification:
[`filterPSI()`](https://nuno-agostinho.github.io/psichomics/reference/filterPSI.md),
[`getSplicingEventTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventTypes.md),
[`listSplicingAnnotations()`](https://nuno-agostinho.github.io/psichomics/reference/listSplicingAnnotations.md),
[`loadAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/loadAnnotation.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md)

## Examples

``` r
# Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
annot <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")

quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
#>                                  Normal 1  Normal 2  Normal 3  Cancer 1
#> SE_1_+_32_35_37_38_ACTN1        0.5000000 0.5000000 0.5000000 0.5000000
#> SE_2_+_32_35_37_38_ALDOA        0.7402597 0.5500000 0.8750000 0.6756757
#> SE_3_+_32_35_37_38_ANXA6        0.4255319 0.5121951 0.3684211 0.4418605
#> MXE_1_+_32_35_37_38_40_42_ACTN1 0.5000000 0.5000000 0.5000000 0.5000000
#> MXE_2_+_32_35_37_38_40_42_ALDOA 0.7500000 0.6000000 0.8803828 0.6500000
#> MXE_3_+_32_35_37_38_40_42_ANXA6 0.2597403 0.4210526 0.1207729 0.3589744
#>                                  Cancer 2  Cancer 3
#> SE_1_+_32_35_37_38_ACTN1        0.5000000 0.5000000
#> SE_2_+_32_35_37_38_ALDOA        0.5744681 0.6037736
#> SE_3_+_32_35_37_38_ANXA6        0.4509804 0.4166667
#> MXE_1_+_32_35_37_38_40_42_ACTN1 0.5000000 0.5000000
#> MXE_2_+_32_35_37_38_40_42_ALDOA 0.5306122 0.5918367
#> MXE_3_+_32_35_37_38_40_42_ANXA6 0.4600000 0.3846154
```
