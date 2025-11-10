# Calculate the contribution of PCA loadings to the selected principal components

Total contribution of a variable is calculated as per
`((Cx * Ex) + (Cy * Ey))/(Ex + Ey)`, where:

- `Cx` and `Cy` are the contributions of a variable to principal
  components `x` and `y`

- `Ex` and `Ey` are the eigenvalues of principal components `x` and `y`

## Usage

``` r
calculateLoadingsContribution(pca, pcX = 1, pcY = 2)
```

## Source

<http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/>

## Arguments

- pca:

  `prcomp` object

- pcX:

  Character: name of the X axis of interest from the PCA

- pcY:

  Character: name of the Y axis of interest from the PCA

## Value

Data frame containing the correlation between variables and selected
principal components and the contribution of variables to the selected
principal components (both individual and total contribution)

## See also

Other functions to analyse principal components:
[`performPCA()`](https://nuno-agostinho.github.io/psichomics/reference/performPCA.md),
[`plotPCA()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCA.md),
[`plotPCAvariance()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCAvariance.md)

## Examples

``` r
pca <- performPCA(USArrests)
calculateLoadingsContribution(pca)
#>          Rank     Gene PC1 loading PC2 loading Contribution to PC1 (%)
#> Assault     1  Assault  0.99522128 -0.05876003              99.0465399
#> UrbanPop    2 UrbanPop  0.04633575  0.97685748               0.2147001
#> Rape        3     Rape  0.07515550  0.20071807               0.5648349
#> Murder      4   Murder  0.04170432 -0.04482166               0.1739250
#>          Contribution to PC2 (%) Contribution to PC1 and PC2 (%)
#> Assault                0.3452741                      96.2825574
#> UrbanPop              95.4250536                       2.8809248
#> Rape                   4.0287742                       0.6618374
#> Murder                 0.2008981                       0.1746804
```
