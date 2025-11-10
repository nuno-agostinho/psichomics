# Create a scatterplot from a PCA object

Create a scatterplot from a PCA object

## Usage

``` r
plotPCA(
  pca,
  pcX = 1,
  pcY = 2,
  groups = NULL,
  individuals = TRUE,
  loadings = FALSE,
  nLoadings = NULL
)
```

## Arguments

- pca:

  `prcomp` object

- pcX:

  Character: name of the X axis of interest from the PCA

- pcY:

  Character: name of the Y axis of interest from the PCA

- groups:

  Matrix: groups to plot indicating the index of interest of the samples
  (use clinical or sample groups)

- individuals:

  Boolean: plot PCA individuals

- loadings:

  Boolean: plot PCA loadings/rotations

- nLoadings:

  Integer: Number of variables to plot, ordered by those that most
  contribute to selected principal components (this allows for faster
  performance as only the most contributing variables are rendered); if
  `NULL`, all variables are plotted

## Value

Scatterplot as an `highchart` object

## See also

Other functions to analyse principal components:
[`calculateLoadingsContribution()`](https://nuno-agostinho.github.io/psichomics/reference/calculateLoadingsContribution.md),
[`performPCA()`](https://nuno-agostinho.github.io/psichomics/reference/performPCA.md),
[`plotPCAvariance()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCAvariance.md)

## Examples

``` r
pca <- prcomp(USArrests, scale=TRUE)
plotPCA(pca)
#> Error: unable to find an inherited method for function ‘plotPCA’ for signature ‘object = "prcomp"’
plotPCA(pca, pcX=2, pcY=3)
#> Error: unable to find an inherited method for function ‘plotPCA’ for signature ‘object = "prcomp"’

# Plot both individuals and loadings
plotPCA(pca, pcX=2, pcY=3, loadings=TRUE)
#> Error: unable to find an inherited method for function ‘plotPCA’ for signature ‘object = "prcomp"’

# Only plot loadings
plotPCA(pca, pcX=2, pcY=3, loadings=TRUE, individuals=FALSE)
#> Error: unable to find an inherited method for function ‘plotPCA’ for signature ‘object = "prcomp"’
```
