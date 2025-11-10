# Basic statistics performed on data

Variance and median of each group. If data has 2 groups, also calculates
the delta variance and delta median.

## Usage

``` r
basicStats(data, groups)
```

## Arguments

- data:

  Numeric, data frame or matrix: gene expression data or alternative
  splicing event quantification values (sample names are based on their
  `names` or `colnames`)

- groups:

  List of sample names or vector containing the group name per `data`
  value (read Details); if `NULL` or a character vector of length 1,
  `data` values are considered from the same group

## Value

HTML elements
