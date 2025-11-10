# Levene's test

Performs a Levene's test to assess the equality of variances

## Usage

``` r
leveneTest(x, g, centers = median)
```

## Arguments

- x:

  Numeric vector or list of numeric vectors: non-numeric elements of a
  list will be coerced with a warning

- g:

  Vector or factor: groups of elements in `x` (ignored with a warning if
  `x` is a list)

- centers:

  Function used to calculate how much values spread; for instance,
  `median` (default) or `mean`

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the test statistic with a name describing it.

- p.value:

  the p-value for the test.

- method:

  the type of test applied.

- data.name:

  a character string giving the names of the data.

## Details

The implementation of this function is based on
`car:::leveneTest.default` with a more standard result.

## Examples

``` r
vals <- sample(30, replace=TRUE)
group <- lapply(list("A", "B", "C"), rep, 10)
group <- unlist(group)
psichomics:::leveneTest(vals, group)
#> 
#>  Levene's test (using the median)
#> 
#> data:  vals and group
#> W = 0.97808, p-value = 0.389
#> 

## Using Levene's test based on the mean
psichomics:::leveneTest(vals, group, mean)
#> 
#>  Levene's test (using the mean)
#> 
#> data:  vals and group
#> W = 1.3031, p-value = 0.2882
#> 
```
