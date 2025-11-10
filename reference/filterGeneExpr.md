# Filter genes based on their expression

Uses [`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
to determine genes with sufficiently large counts to retain for
statistical analysis.

## Usage

``` r
filterGeneExpr(
  geneExpr,
  minMean = 0,
  maxMean = Inf,
  minVar = 0,
  maxVar = Inf,
  minCounts = 10,
  minTotalCounts = 15
)
```

## Arguments

- geneExpr:

  Data frame or matrix: gene expression

- minMean:

  Numeric: minimum of read count mean per gene

- maxMean:

  Numeric: maximum of read count mean per gene

- minVar:

  Numeric: minimum of read count variance per gene

- maxVar:

  Numeric: maximum of read count variance per gene

- minCounts:

  Numeric: minimum number of read counts per gene for a worthwhile
  number of samples (check
  [`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html) for
  more information)

- minTotalCounts:

  Numeric: minimum total number of read counts per gene

## Value

Boolean vector indicating which genes have sufficiently large counts

## See also

Other functions for gene expression pre-processing:
[`convertGeneIdentifiers()`](https://nuno-agostinho.github.io/psichomics/reference/convertGeneIdentifiers.md),
[`normaliseGeneExpression()`](https://nuno-agostinho.github.io/psichomics/reference/normaliseGeneExpression.md),
[`plotGeneExprPerSample()`](https://nuno-agostinho.github.io/psichomics/reference/plotGeneExprPerSample.md),
[`plotLibrarySize()`](https://nuno-agostinho.github.io/psichomics/reference/plotLibrarySize.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md)

## Examples

``` r
geneExpr <- readFile("ex_gene_expression.RDS")

# Add some genes with low expression
geneExpr <- rbind(geneExpr,
                  lowReadGene1=c(rep(4:5, 10)),
                  lowReadGene2=c(rep(5:1, 10)),
                  lowReadGene3=c(rep(10:1, 10)),
                  lowReadGene4=c(rep(7:8, 10)))
#> Warning: number of columns of result, 10, is not a multiple of vector length 20 of arg 2
#> Warning: number of columns of result, 10, is not a multiple of vector length 50 of arg 3
#> Warning: number of columns of result, 10, is not a multiple of vector length 100 of arg 4
#> Warning: number of columns of result, 10, is not a multiple of vector length 20 of arg 5

# Filter out genes with low reads across samples
geneExpr[filterGeneExpr(geneExpr), ]
#>                Normal 1 Normal 2 Normal 3 Normal 4 Normal 5 Cancer 1 Cancer 2
#> ACTN1|87       15810.00 10454.00  6247.00 23657.00 17107.00 13509.00  7038.00
#> ADAM15|8751    18894.97  6483.95  5463.78  9740.44  5042.07  4866.82  4901.09
#> AKAP8L|26993    3211.00  2661.00  1322.00  3322.00  5009.00  2366.00  1358.00
#> AKR1A1|10327    6413.00  3170.00  2018.00  8364.00  5797.00  4687.00  3901.00
#> ALDOA|226      69936.00 97762.00 26714.00 72208.00 47405.00 31632.00 34473.00
#> ANXA6|309      20259.00 10964.00  7139.00 16940.00 12995.00 10893.00  8817.00
#> ARAP1|116985   11774.26  4539.73  3070.90  6822.52  7186.83  6104.56  4634.95
#> ATP5C1|509      4869.00  4210.00  4822.00  5879.00  6267.00  4283.00  3530.00
#> BOLA2|552900    2978.44  3153.36  1790.55  1806.70  2535.48  2035.52  1330.20
#> C16orf13|84326  2453.00  1758.00  1144.00  2028.00  2381.00  1355.00   952.00
#>                Cancer 3 Cancer 4 Cancer 5
#> ACTN1|87        9426.00 21930.00 27954.00
#> ADAM15|8751     7038.25  3867.00  9977.15
#> AKAP8L|26993    2958.00  1946.00  1602.00
#> AKR1A1|10327    3351.00  3670.00  3807.00
#> ALDOA|226      28357.00 22264.00 53072.00
#> ANXA6|309      11703.00  5670.00 23697.00
#> ARAP1|116985    6058.86  4088.99  6724.75
#> ATP5C1|509      7696.00  5527.00  8987.00
#> BOLA2|552900    2615.29  1049.42  2471.20
#> C16orf13|84326  1214.00   904.00  1703.00
```
