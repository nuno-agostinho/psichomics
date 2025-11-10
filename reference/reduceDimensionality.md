# Reduce dimensionality after processing missing values from data frame

Reduce dimensionality after processing missing values from data frame

## Usage

``` r
reduceDimensionality(
  data,
  type = c("pca", "ica"),
  center = TRUE,
  scale. = FALSE,
  naTolerance = NULL,
  missingValues = round(0.05 * ncol(data)),
  ...
)
```

## Arguments

- data:

  Data frame: data

- type:

  Character: dimensionality reduction technique (`pca` or `ica`)

- center:

  either a logical value or numeric-alike vector of length equal to the
  number of columns of `x`, where ‘numeric-alike’ means that
  [`as.numeric`](https://rdrr.io/r/base/numeric.html)`(.)` will be
  applied successfully if
  [`is.numeric`](https://rdrr.io/r/base/numeric.html)`(.)` is not true.

- scale.:

  Boolean: scale variables?

- naTolerance:

  Integer: percentage of tolerated missing values per column
  (deprecated)

- missingValues:

  Integer: number of tolerated missing values per column to be replaced
  with the mean of the values of that same column

- ...:

  Extra parameters passed to FUN

## Value

PCA result in a `prcomp` object or ICA result object
