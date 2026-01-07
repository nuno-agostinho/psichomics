# Sum columns using an `EList-class` object

Sum columns using an `EList-class` object

## Usage

``` r
# S4 method for class 'EList'
colSums(x, na.rm = FALSE, dims = 1)
```

## Arguments

- x:

  an array of two or more dimensions, containing numeric, complex,
  integer or logical values, or a numeric data frame. For
  [`.colSums()`](https://rdrr.io/r/base/colSums.html) etc, a numeric,
  integer or logical matrix (or vector of length `m * n`).

- na.rm:

  logical. Should missing values (including `NaN`) be omitted from the
  calculations?

- dims:

  integer number: Which dimensions are regarded as ‘rows’ or ‘columns’
  to sum over. For `row*`, the sum or mean is over dimensions
  `dims+1, ...`; for `col*` it is over dimensions `1:dims`.

## Value

Numeric vector with the sum of the columns
