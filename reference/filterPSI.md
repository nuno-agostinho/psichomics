# Filter alternative splicing quantification

Filter alternative splicing quantification

## Usage

``` r
filterPSI(
  psi,
  eventType = NULL,
  eventSubtype = NULL,
  minPSI = -Inf,
  maxPSI = Inf,
  minMedian = -Inf,
  maxMedian = Inf,
  minLogVar = -Inf,
  maxLogVar = Inf,
  minRange = -Inf,
  maxRange = Inf
)
```

## Arguments

- psi:

  Data frame or matrix: alternative splicing quantification

- eventType:

  Character: filter data based on event type; check all event types
  available by using `getSplicingEventTypes(psi)`, where `psi` is the
  alternative splicing quantification data; if `eventType = NULL`,
  events are not filtered by event type

- eventSubtype:

  Character: filter data based on event subtype; check all event
  subtypes available in your data by using
  `unique(getSplicingEventData(psi)$subtype)`, where `psi` is the
  alternative splicing quantification data; if `eventSubtype = NULL`,
  events are not filtered by event subtype

- minPSI:

  Numeric: minimum PSI value

- maxPSI:

  Numeric: maximum PSI value

- minMedian:

  Numeric: minimum median PSI per splicing event

- maxMedian:

  Numeric: maximum median PSI per splicing event

- minLogVar:

  Numeric: minimum log10(PSI variance) per splicing event

- maxLogVar:

  Numeric: maximum log10(PSI variance) per splicing event

- minRange:

  Numeric: minimum PSI range across samples per splicing event

- maxRange:

  Numeric: maximum PSI range across samples per splicing event

## Value

Boolean vector indicating which splicing events pass the thresholds

## See also

Other functions for PSI quantification:
[`getSplicingEventTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventTypes.md),
[`listSplicingAnnotations()`](https://nuno-agostinho.github.io/psichomics/reference/listSplicingAnnotations.md),
[`loadAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/loadAnnotation.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md),
[`quantifySplicing()`](https://nuno-agostinho.github.io/psichomics/reference/quantifySplicing.md)

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
# Filter PSI
psi[filterPSI(psi, minMedian=0.05, maxMedian=0.95, minRange=0.15), ]
#>                                  Normal 1  Normal 2  Normal 3  Cancer 1
#> SE_2_+_32_35_37_38_ALDOA        0.7402597 0.5500000 0.8750000 0.6756757
#> MXE_2_+_32_35_37_38_40_42_ALDOA 0.7500000 0.6000000 0.8803828 0.6500000
#> MXE_3_+_32_35_37_38_40_42_ANXA6 0.2597403 0.4210526 0.1207729 0.3589744
#>                                  Cancer 2  Cancer 3
#> SE_2_+_32_35_37_38_ALDOA        0.5744681 0.6037736
#> MXE_2_+_32_35_37_38_40_42_ALDOA 0.5306122 0.5918367
#> MXE_3_+_32_35_37_38_40_42_ANXA6 0.4600000 0.3846154
```
