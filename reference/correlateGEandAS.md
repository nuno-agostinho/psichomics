# Correlate gene expression data against alternative splicing quantification

Test for association between paired samples' gene expression (for any
genes of interest) and alternative splicing quantification.

## Usage

``` r
correlateGEandAS(geneExpr, psi, gene, ASevents = NULL, ...)
```

## Arguments

- geneExpr:

  Matrix or data frame: gene expression data

- psi:

  Matrix or data frame: alternative splicing quantification data

- gene:

  Character: gene symbol for genes of interest

- ASevents:

  Character: alternative splicing events to correlate with gene
  expression of a gene (if `NULL`, the events will be automatically
  retrieved from the given gene)

- ...:

  Extra parameters passed to
  [cor.test](https://rdrr.io/r/stats/cor.test.html)

## Value

List of correlations where each element contains:

- `eventID`:

  Alternative splicing event identifier

- `cor`:

  Correlation between gene expression and alternative splicing
  quantification of one alternative splicing event

- `geneExpr`:

  Gene expression for the selected gene

- `psi`:

  Alternative splicing quantification for the alternative splicing event

## See also

Other functions to correlate gene expression and alternative splicing:
[`[.GEandAScorrelation()`](https://nuno-agostinho.github.io/psichomics/reference/plot.GEandAScorrelation.md)

## Examples

``` r
annot <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")
psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 

geneExpr <- readFile("ex_gene_expression.RDS")
correlateGEandAS(geneExpr, psi, "ALDOA")
#> ================================================================================
#> SE_2_+_32_35_37_38_ALDOA splicing event
#> ALDOA|226 gene expression
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  exprNum and eventPSInum
#> t = -0.7542, df = 4, p-value = 0.4927
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  -0.9051981  0.6427793
#> sample estimates:
#>        cor 
#> -0.3528456 
#> 
#> ================================================================================
#> MXE_2_+_32_35_37_38_40_42_ALDOA splicing event
#> ALDOA|226 gene expression
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  exprNum and eventPSInum
#> t = -0.26642, df = 4, p-value = 0.8031
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  -0.8522745  0.7610748
#> sample estimates:
#>        cor 
#> -0.1320457 
#> 
```
