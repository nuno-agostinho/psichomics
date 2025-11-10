# Filter and normalise gene expression

Gene expression is filtered and normalised in the following steps:

- Filter gene expression;

- Normalise gene expression with
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html);

- If `performVoom = FALSE`, compute counts per million (CPM) using
  [`cpm`](https://rdrr.io/pkg/edgeR/man/cpm.html) and log2-transform
  values if `log2transform = TRUE`;

- If `performVoom = TRUE`, use
  [`voom`](https://rdrr.io/pkg/limma/man/voom.html) to compute log2-CPM,
  quantile-normalise (if `method = "quantile"`) and estimate
  mean-variance relationship to calculate observation-level weights.

## Usage

``` r
normaliseGeneExpression(
  geneExpr,
  geneFilter = NULL,
  method = "TMM",
  p = 0.75,
  log2transform = TRUE,
  priorCount = 0.25,
  performVoom = FALSE
)

normalizeGeneExpression(
  geneExpr,
  geneFilter = NULL,
  method = "TMM",
  p = 0.75,
  log2transform = TRUE,
  priorCount = 0.25,
  performVoom = FALSE
)
```

## Arguments

- geneExpr:

  Matrix or data frame: gene expression

- geneFilter:

  Boolean: filtered genes (if `NULL`, skip filtering)

- method:

  Character: normalisation method, including `TMM`, `RLE`,
  `upperquartile`, `none` or `quantile` (see Details)

- p:

  numeric value between 0 and 1 specifying which quantile of the counts
  should be used by `method="upperquartile"`.

- log2transform:

  Boolean: perform log2-transformation?

- priorCount:

  Average count to add to each observation to avoid zeroes after
  log-transformation

- performVoom:

  Boolean: perform mean-variance modelling (using
  [`voom`](https://rdrr.io/pkg/limma/man/voom.html))?

## Value

Filtered and normalised gene expression

## Details

[`edgeR::calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
will be used to normalise gene expression if `method` is `TMM`, `RLE`,
`upperquartile` or `none`. If `performVoom = TRUE`,
[`voom`](https://rdrr.io/pkg/limma/man/voom.html) will only normalise if
`method = "quantile"`.

Available normalisation methods:

- `TMM` is recommended for most RNA-seq data where more than half of the
  genes are believed not differentially expressed between any pair of
  samples;

- `RLE` calculates the median library from the geometric mean of all
  columns and the median ratio of each sample to the median library is
  taken as the scale factor;

- `upperquartile` calculates the scale factors from a given quantile of
  the counts for each library, after removing genes with zero counts in
  all libraries;

- `quantile` forces the entire empirical distribution of each column to
  be identical (only performed if `performVoom = TRUE`).

## See also

Other functions for gene expression pre-processing:
[`convertGeneIdentifiers()`](https://nuno-agostinho.github.io/psichomics/reference/convertGeneIdentifiers.md),
[`filterGeneExpr()`](https://nuno-agostinho.github.io/psichomics/reference/filterGeneExpr.md),
[`plotGeneExprPerSample()`](https://nuno-agostinho.github.io/psichomics/reference/plotGeneExprPerSample.md),
[`plotLibrarySize()`](https://nuno-agostinho.github.io/psichomics/reference/plotLibrarySize.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md)

## Examples

``` r
geneExpr <- readFile("ex_gene_expression.RDS")
normaliseGeneExpression(geneExpr)
#>                Normal 1 Normal 2 Normal 3 Normal 4 Normal 5 Cancer 1 Cancer 2
#> ACTN1|87       16.44383 16.78304 16.48223 17.24207 17.07503 17.10327 16.49269
#> ADAM15|8751    16.70099 16.09395 16.28898 15.96189 15.31259 15.63044 15.97064
#> AKAP8L|26993   14.14424 14.80911 14.24193 14.41006 15.30310 14.58999 14.11916
#> AKR1A1|10327   15.14212 15.06161 14.85207 15.74210 15.51387 15.57613 15.64140
#> ALDOA|226      18.58900 20.00823 18.57857 18.85194 18.54547 18.33072 18.78489
#> ANXA6|309      16.80155 16.85176 16.67479 16.76024 16.67841 16.79276 16.81781
#> ARAP1|116985   16.01864 15.57970 15.45777 15.44823 15.82390 15.95734 15.89010
#> ATP5C1|509     14.74478 15.47092 16.10871 15.23351 15.62633 15.44610 15.49723
#> BOLA2|552900   14.03579 15.05401 14.67957 13.53149 14.32091 14.37295 14.08933
#> C16orf13|84326 13.75582 14.21114 14.03332 13.69816 14.23023 13.78593 13.60679
#>                Cancer 3 Cancer 4 Cancer 5
#> ACTN1|87       16.81314 18.48316 17.61675
#> ADAM15|8751    16.39172 15.97958 16.13043
#> AKAP8L|26993   15.14118 14.98893 13.49192
#> AKR1A1|10327   15.32114 15.90415 14.74053
#> ALDOA|226      18.40211 18.50496 18.54164
#> ANXA6|309      17.12530 16.53170 17.37841
#> ARAP1|116985   16.17556 16.06011 15.56131
#> ATP5C1|509     16.52061 16.49485 15.97965
#> BOLA2|552900   14.96354 14.09810 14.11716
#> C16orf13|84326 13.85646 13.88293 13.58011
```
