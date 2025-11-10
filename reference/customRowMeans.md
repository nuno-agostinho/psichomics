# Calculate statistics for each row or column of a matrix

Calculate statistics for each row or column of a matrix

## Usage

``` r
customRowMeans(mat, na.rm = FALSE, fast = FALSE)

customRowMedians(mat, na.rm = FALSE, fast = FALSE)

customRowVars(mat, na.rm = FALSE, fast = FALSE)

customRowMins(mat, na.rm = FALSE, fast = FALSE)

customRowMaxs(mat, na.rm = FALSE, fast = FALSE)

customRowRanges(mat, na.rm = FALSE, fast = FALSE)

customColMedians(mat, na.rm = FALSE, fast = FALSE)
```

## Arguments

- mat:

  Matrix

- na.rm:

  Boolean: remove missing values (`NA`)?

- fast:

  Boolean: use `Rfast` functions? They may return different results from
  R built-in functions

## Value

Vector of selected statistic

## Examples

``` r
df <- rbind("Gene 1"=c(3, 5, 7), "Gene 2"=c(8, 2, 4), "Gene 3"=c(9:11))
psichomics:::customRowMeans(df)
#>    Gene 1    Gene 2    Gene 3 
#>  5.000000  4.666667 10.000000 
psichomics:::customRowVars(df, fast=TRUE)
#>   Gene 1   Gene 2   Gene 3 
#> 4.000000 9.333333 1.000000 
```
