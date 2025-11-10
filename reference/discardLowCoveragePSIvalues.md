# Remove alternative splicing quantification values based on coverage

Remove alternative splicing quantification values based on coverage

## Usage

``` r
discardLowCoveragePSIvalues(
  psi,
  minReads = 10,
  vasttoolsScoresToDiscard = c("VLOW", "N")
)
```

## Arguments

- psi:

  Data frame or matrix: alternative splicing quantification

- minReads:

  Currently this argument does nothing

- vasttoolsScoresToDiscard:

  Character: if you are using inclusion levels from VAST-TOOLS, filter
  the data based on quality scores for read coverage, e.g. use
  `vasttoolsScoresToDiscard = c("SOK", "OK", "LOW")` to only keep events
  with good read coverage (by default, events are not filtered based on
  quality scores); read <https://github.com/vastgroup/vast-tools> for
  more information on VAST-TOOLS quality scores

## Value

Alternative splicing quantification data with missing values for any
values with insufficient coverage
