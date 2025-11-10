# Filters the events with valid elements according to the given validator

Filters the events with valid elements according to the given validator

## Usage

``` r
getValidEvents(event, validator, areMultipleExonsValid = FALSE)
```

## Arguments

- event:

  Data.frame containing only one event with at least 7 columns as
  retrieved from the alternative splicing annotation files from MISO
  (GFF3 files)

- validator:

  Character: valid elements for each event

- areMultipleExonsValid:

  Boolean: consider runs of exons as valid when comparing with the
  validator? Default is `FALSE` (see details)

## Value

Data.frame with valid events

## Details

`areMultipleExonsValid` allows to consider runs of exons (i.e. sequences
where `exon` occurs consecutively) as valid when comparing based on the
`validator`. For example, if `validator = c("gene", "mRNA", "exon")` and
`areMultipleExonsValid = FALSE`, the event
`c("gene", "mRNA", "exon", "exon")` is not valid as it has one
additional exon. If `areMultipleExonsValid = TRUE`, the same event would
be valid.

## Examples

``` r
event <- read.table(text = "
 chr1 SE gene 17233 18061  .  -  .
 chr1 SE dkfd 00000 30000  .  -  .
 chr1 SE mRNA 17233 18061  .  -  .
 chr1 SE exon 17233 17368  .  -  .
 chr1 SE exon 17526 17742  .  -  .
 chr1 SE exon 17915 18061  .  -  .
 chr1 SE mRNA 17233 18061  .  -  .
 chr1 SE exon 17233 17368  .  -  .
 chr1 SE exon 17915 18061  .  -  .
 chr1 SE gene 17233 18061  .  -  .
 chr1 SE mRNA 17233 18061  .  -  .
 chr1 SE exon 17233 17368  .  -  .
 chr1 SE exon 17606 17742  .  -  .
 chr1 SE exon 17915 18061  .  -  .
 chr1 SE mRNA 17233 18061  .  -  .
 chr1 SE exon 17233 17368  .  -  .
 chr1 SE exon 17915 18061  .  -  .
")
validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 2))
psichomics:::getValidEvents(event, validator)
#>      V1 V2   V3    V4    V5 V6 V7 V8
#> 10 chr1 SE gene 17233 18061  .  -  .
#> 11 chr1 SE mRNA 17233 18061  .  -  .
#> 12 chr1 SE exon 17233 17368  .  -  .
#> 13 chr1 SE exon 17606 17742  .  -  .
#> 14 chr1 SE exon 17915 18061  .  -  .
#> 15 chr1 SE mRNA 17233 18061  .  -  .
#> 16 chr1 SE exon 17233 17368  .  -  .
#> 17 chr1 SE exon 17915 18061  .  -  .
```
