# Perform statistical analysis on a given splicing event

Perform statistical analyses on a given vector containing elements from
different groups

## Usage

``` r
singleDiffAnalyses(
  vector,
  group,
  threshold = 1,
  step = 100,
  analyses = c("wilcoxRankSum", "ttest", "kruskal", "levene", "fligner")
)
```

## Arguments

- vector:

  Numeric

- group:

  Character: group of each element in the vector

- threshold:

  Integer: minimum number of values per group

- step:

  Numeric: number of events before the progress bar is updated (a bigger
  number allows for a faster execution)

- analyses:

  Character: analyses to perform (see Details)

## Value

A row from a data frame with the results

## Details

The following statistical analyses may be performed by including the
respective string in the `analysis` argument:

- `ttest` - Unpaired t-test (2 groups)

- `wilcoxRankSum` - Wilcoxon Rank Sum test (2 groups)

- `kruskal` - Kruskal test (2 or more groups)

- `levene` - Levene's test (2 or more groups)

- `fligner` - Fligner-Killeen test (2 or more groups)
