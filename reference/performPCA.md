# Perform principal component analysis after processing missing values

Perform principal component analysis after processing missing values

## Usage

``` r
performPCA(
  data,
  center = TRUE,
  scale. = FALSE,
  missingValues = round(0.05 * nrow(data)),
  ...
)
```

## Arguments

- data:

  an optional data frame (or similar: see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula `formula`. By default the variables are
  taken from `environment(formula)`.

- center:

  a logical value indicating whether the variables should be shifted to
  be zero centered. Alternately, a vector of length equal the number of
  columns of `x` can be supplied. The value is passed to `scale`.

- scale.:

  a logical value indicating whether the variables should be scaled to
  have unit variance before the analysis takes place. The default is
  `FALSE` for consistency with S, but in general scaling is advisable.
  Alternatively, a vector of length equal the number of columns of `x`
  can be supplied. The value is passed to
  [`scale`](https://rdrr.io/r/base/scale.html).

- missingValues:

  Integer: number of tolerated missing values per column to be replaced
  with the mean of the values of that same column

- ...:

  Arguments passed on to
  [`stats::prcomp`](https://rdrr.io/r/stats/prcomp.html)

## Value

PCA result in a `prcomp` object

## See also

Other functions to analyse principal components:
[`calculateLoadingsContribution()`](https://nuno-agostinho.github.io/psichomics/reference/calculateLoadingsContribution.md),
[`plotPCA()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCA.md),
[`plotPCAvariance()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCAvariance.md)

## Examples

``` r
performPCA(USArrests)
#> Standard deviations (1, .., p=4):
#> [1] 83.732400 14.212402  6.489426  2.482790
#> 
#> Rotation (n x k) = (4 x 4):
#>                 PC1         PC2         PC3         PC4
#> Murder   0.04170432 -0.04482166  0.07989066 -0.99492173
#> Assault  0.99522128 -0.05876003 -0.06756974  0.03893830
#> UrbanPop 0.04633575  0.97685748 -0.20054629 -0.05816914
#> Rape     0.07515550  0.20071807  0.97408059  0.07232502
```
