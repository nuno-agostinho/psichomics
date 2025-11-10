# Perform and display statistical analysis

Includes interface containing the results

## Usage

``` r
wilcox(data, groups, stat = NULL)

ttest(data, groups, stat = NULL)

levene(data, groups, stat = NULL)

fligner(data, groups, stat = NULL)

kruskal(data, groups, stat = NULL)

fisher(data, groups)

spearman(data, groups)
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

- stat:

  Data frame or matrix: values of the analyses to be performed (if
  `NULL`, the analyses will be performed)

## Value

HTML elements

## Details

- `ttest`: unpaired t-test

- `wilcox`: Wilcoxon test

- `levene`: Levene's test

- `fligner`: Fligner-Killeen test

- `kruskal`: Kruskal test

- `fisher`: Fisher's exact test

- `spearman`: Spearman's test
