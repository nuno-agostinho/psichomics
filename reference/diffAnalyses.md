# Perform statistical analyses

Perform statistical analyses

## Usage

``` r
diffAnalyses(
  data,
  groups = NULL,
  analyses = c("wilcoxRankSum", "ttest", "kruskal", "levene", "fligner"),
  pvalueAdjust = "BH",
  geneExpr = NULL,
  inputID = "sparklineInput"
)
```

## Arguments

- data:

  Data frame or matrix: gene expression or alternative splicing
  quantification

- groups:

  Named list of characters (containing elements belonging to each group)
  or character vector (containing the group of each individual sample);
  if `NULL`, sample types are used instead when available, e.g. normal,
  tumour and metastasis

- analyses:

  Character: statistical tests to perform (see Details)

- pvalueAdjust:

  Character: method used to adjust p-values (see Details)

- geneExpr:

  Character: name of the gene expression dataset (only required for
  density sparklines available in the interactive mode)

- inputID:

  Character: identifier of input to get attributes of clicked event
  (Shiny only)

## Value

Table of statistical analyses

## Details

The following statistical analyses may be performed simultaneously via
the `analysis` argument:

- `ttest` - Unpaired t-test (2 groups)

- `wilcoxRankSum` - Wilcoxon Rank Sum test (2 groups)

- `kruskal` - Kruskal test (2 or more groups)

- `levene` - Levene's test (2 or more groups)

- `fligner` - Fligner-Killeen test (2 or more groups)

- `density` - Sample distribution per group (only usable through the
  visual interface)

The following p-value adjustment methods are supported via the
`pvalueAdjust` argument:

- `none`: do not adjust p-values

- `BH`: Benjamini-Hochberg's method (false discovery rate)

- `BY`: Benjamini-Yekutieli's method (false discovery rate)

- `bonferroni`: Bonferroni correction (family-wise error rate)

- `holm`: Holm's method (family-wise error rate)

- `hochberg`: Hochberg's method (family-wise error rate)

- `hommel`: Hommel's method (family-wise error rate)

## See also

Other functions to perform and plot differential analyses:
[`plotDistribution()`](https://nuno-agostinho.github.io/psichomics/reference/plotDistribution.md)

## Examples

``` r
# Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
eventType <- c("SE", "MXE")
annot <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")

psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
group <- c(rep("Normal", 3), rep("Tumour", 3))
diffAnalyses(psi, group)
#> 
#> Time difference of 0.0225 secs
#> 
#> Time difference of 0.00589 secs
#> 
#> Time difference of 0.000306 secs
#> 
#> Time difference of 0.000975 secs
#>                                                    Event type Chromosome Strand
#> SE_1_+_32_35_37_38_ACTN1                    Skipped exon (SE)          1      +
#> MXE_1_+_32_35_37_38_40_42_ACTN1 Mutually exclusive exon (MXE)          1      +
#> SE_2_+_32_35_37_38_ALDOA                    Skipped exon (SE)          2      +
#> SE_3_+_32_35_37_38_ANXA6                    Skipped exon (SE)          3      +
#> MXE_2_+_32_35_37_38_40_42_ALDOA Mutually exclusive exon (MXE)          2      +
#> MXE_3_+_32_35_37_38_40_42_ANXA6 Mutually exclusive exon (MXE)          3      +
#>                                  Gene Survival by PSI cutoff Optimal PSI cutoff
#> SE_1_+_32_35_37_38_ACTN1        ACTN1                   <NA>               <NA>
#> MXE_1_+_32_35_37_38_40_42_ACTN1 ACTN1                   <NA>               <NA>
#> SE_2_+_32_35_37_38_ALDOA        ALDOA                   <NA>               <NA>
#> SE_3_+_32_35_37_38_ANXA6        ANXA6                   <NA>               <NA>
#> MXE_2_+_32_35_37_38_40_42_ALDOA ALDOA                   <NA>               <NA>
#> MXE_3_+_32_35_37_38_40_42_ANXA6 ANXA6                   <NA>               <NA>
#>                                 Log-rank p-value Samples (Normal)
#> SE_1_+_32_35_37_38_ACTN1                    <NA>                3
#> MXE_1_+_32_35_37_38_40_42_ACTN1             <NA>                3
#> SE_2_+_32_35_37_38_ALDOA                    <NA>                3
#> SE_3_+_32_35_37_38_ANXA6                    <NA>                3
#> MXE_2_+_32_35_37_38_40_42_ALDOA             <NA>                3
#> MXE_3_+_32_35_37_38_40_42_ANXA6             <NA>                3
#>                                 Samples (Tumour) Wilcoxon statistic
#> SE_1_+_32_35_37_38_ACTN1                       3                4.5
#> MXE_1_+_32_35_37_38_40_42_ACTN1                3                4.5
#> SE_2_+_32_35_37_38_ALDOA                       3                6.0
#> SE_3_+_32_35_37_38_ANXA6                       3                4.0
#> MXE_2_+_32_35_37_38_40_42_ALDOA                3                8.0
#> MXE_3_+_32_35_37_38_40_42_ANXA6                3                2.0
#>                                 Wilcoxon p-value Wilcoxon p-value (BH adjusted)
#> SE_1_+_32_35_37_38_ACTN1                     NaN                            NaN
#> MXE_1_+_32_35_37_38_40_42_ACTN1              NaN                            NaN
#> SE_2_+_32_35_37_38_ALDOA                     0.7                      0.9333333
#> SE_3_+_32_35_37_38_ANXA6                     1.0                      1.0000000
#> MXE_2_+_32_35_37_38_40_42_ALDOA              0.2                      0.8000000
#> MXE_3_+_32_35_37_38_40_42_ANXA6              0.4                      0.8000000
#>                                 Wilcoxon null value Wilcoxon alternative
#> SE_1_+_32_35_37_38_ACTN1                          0            two.sided
#> MXE_1_+_32_35_37_38_40_42_ACTN1                   0            two.sided
#> SE_2_+_32_35_37_38_ALDOA                          0            two.sided
#> SE_3_+_32_35_37_38_ANXA6                          0            two.sided
#> MXE_2_+_32_35_37_38_40_42_ALDOA                   0            two.sided
#> MXE_3_+_32_35_37_38_40_42_ANXA6                   0            two.sided
#>                                                                   Wilcoxon method
#> SE_1_+_32_35_37_38_ACTN1        Wilcoxon rank sum test with continuity correction
#> MXE_1_+_32_35_37_38_40_42_ACTN1 Wilcoxon rank sum test with continuity correction
#> SE_2_+_32_35_37_38_ALDOA                             Wilcoxon rank sum exact test
#> SE_3_+_32_35_37_38_ANXA6                             Wilcoxon rank sum exact test
#> MXE_2_+_32_35_37_38_40_42_ALDOA                      Wilcoxon rank sum exact test
#> MXE_3_+_32_35_37_38_40_42_ANXA6                      Wilcoxon rank sum exact test
#>                                                   Wilcoxon data name
#> SE_1_+_32_35_37_38_ACTN1        vector[typeOne] and vector[!typeOne]
#> MXE_1_+_32_35_37_38_40_42_ACTN1 vector[typeOne] and vector[!typeOne]
#> SE_2_+_32_35_37_38_ALDOA        vector[typeOne] and vector[!typeOne]
#> SE_3_+_32_35_37_38_ANXA6        vector[typeOne] and vector[!typeOne]
#> MXE_2_+_32_35_37_38_40_42_ALDOA vector[typeOne] and vector[!typeOne]
#> MXE_3_+_32_35_37_38_40_42_ANXA6 vector[typeOne] and vector[!typeOne]
#>                                 Kruskal statistic Kruskal parameter
#> SE_1_+_32_35_37_38_ACTN1                      NaN                 1
#> MXE_1_+_32_35_37_38_40_42_ACTN1               NaN                 1
#> SE_2_+_32_35_37_38_ALDOA               0.42857143                 1
#> SE_3_+_32_35_37_38_ANXA6               0.04761905                 1
#> MXE_2_+_32_35_37_38_40_42_ALDOA        2.33333333                 1
#> MXE_3_+_32_35_37_38_40_42_ANXA6        1.19047619                 1
#>                                 Kruskal p-value Kruskal p-value (BH adjusted)
#> SE_1_+_32_35_37_38_ACTN1                    NaN                           NaN
#> MXE_1_+_32_35_37_38_40_42_ACTN1             NaN                           NaN
#> SE_2_+_32_35_37_38_ALDOA              0.5126908                     0.6835877
#> SE_3_+_32_35_37_38_ANXA6              0.8272593                     0.8272593
#> MXE_2_+_32_35_37_38_40_42_ALDOA       0.1266305                     0.5065218
#> MXE_3_+_32_35_37_38_40_42_ANXA6       0.2752335                     0.5504670
#>                                               Kruskal method Kruskal data name
#> SE_1_+_32_35_37_38_ACTN1        Kruskal-Wallis rank sum test  vector and group
#> MXE_1_+_32_35_37_38_40_42_ACTN1 Kruskal-Wallis rank sum test  vector and group
#> SE_2_+_32_35_37_38_ALDOA        Kruskal-Wallis rank sum test  vector and group
#> SE_3_+_32_35_37_38_ANXA6        Kruskal-Wallis rank sum test  vector and group
#> MXE_2_+_32_35_37_38_40_42_ALDOA Kruskal-Wallis rank sum test  vector and group
#> MXE_3_+_32_35_37_38_40_42_ANXA6 Kruskal-Wallis rank sum test  vector and group
#>                                 Levene statistic Levene p-value
#> SE_1_+_32_35_37_38_ACTN1                     NaN            NaN
#> MXE_1_+_32_35_37_38_40_42_ACTN1              NaN            NaN
#> SE_2_+_32_35_37_38_ALDOA                1.534423      0.2831735
#> SE_3_+_32_35_37_38_ANXA6                1.898555      0.2403023
#> MXE_2_+_32_35_37_38_40_42_ALDOA         1.102360      0.3529948
#> MXE_3_+_32_35_37_38_40_42_ANXA6         1.453025      0.2944743
#>                                 Levene p-value (BH adjusted) Levene data name
#> SE_1_+_32_35_37_38_ACTN1                                 NaN vector and group
#> MXE_1_+_32_35_37_38_40_42_ACTN1                          NaN vector and group
#> SE_2_+_32_35_37_38_ALDOA                           0.3529948 vector and group
#> SE_3_+_32_35_37_38_ANXA6                           0.3529948 vector and group
#> MXE_2_+_32_35_37_38_40_42_ALDOA                    0.3529948 vector and group
#> MXE_3_+_32_35_37_38_40_42_ANXA6                    0.3529948 vector and group
#>                                                    Levene method
#> SE_1_+_32_35_37_38_ACTN1        Levene's test (using the median)
#> MXE_1_+_32_35_37_38_40_42_ACTN1 Levene's test (using the median)
#> SE_2_+_32_35_37_38_ALDOA        Levene's test (using the median)
#> SE_3_+_32_35_37_38_ANXA6        Levene's test (using the median)
#> MXE_2_+_32_35_37_38_40_42_ALDOA Levene's test (using the median)
#> MXE_3_+_32_35_37_38_40_42_ANXA6 Levene's test (using the median)
#>                                 Fligner-Killeen statistic
#> SE_1_+_32_35_37_38_ACTN1                              NaN
#> MXE_1_+_32_35_37_38_40_42_ACTN1                       NaN
#> SE_2_+_32_35_37_38_ALDOA                         1.042057
#> SE_3_+_32_35_37_38_ANXA6                         1.042057
#> MXE_2_+_32_35_37_38_40_42_ALDOA                  1.042057
#> MXE_3_+_32_35_37_38_40_42_ANXA6                  1.042057
#>                                 Fligner-Killeen parameter
#> SE_1_+_32_35_37_38_ACTN1                                1
#> MXE_1_+_32_35_37_38_40_42_ACTN1                         1
#> SE_2_+_32_35_37_38_ALDOA                                1
#> SE_3_+_32_35_37_38_ANXA6                                1
#> MXE_2_+_32_35_37_38_40_42_ALDOA                         1
#> MXE_3_+_32_35_37_38_40_42_ANXA6                         1
#>                                 Fligner-Killeen p-value
#> SE_1_+_32_35_37_38_ACTN1                            NaN
#> MXE_1_+_32_35_37_38_40_42_ACTN1                     NaN
#> SE_2_+_32_35_37_38_ALDOA                      0.3073435
#> SE_3_+_32_35_37_38_ANXA6                      0.3073435
#> MXE_2_+_32_35_37_38_40_42_ALDOA               0.3073435
#> MXE_3_+_32_35_37_38_40_42_ANXA6               0.3073435
#>                                 Fligner-Killeen p-value (BH adjusted)
#> SE_1_+_32_35_37_38_ACTN1                                          NaN
#> MXE_1_+_32_35_37_38_40_42_ACTN1                                   NaN
#> SE_2_+_32_35_37_38_ALDOA                                    0.3073435
#> SE_3_+_32_35_37_38_ANXA6                                    0.3073435
#> MXE_2_+_32_35_37_38_40_42_ALDOA                             0.3073435
#> MXE_3_+_32_35_37_38_40_42_ANXA6                             0.3073435
#>                                                           Fligner-Killeen method
#> SE_1_+_32_35_37_38_ACTN1        Fligner-Killeen test of homogeneity of variances
#> MXE_1_+_32_35_37_38_40_42_ACTN1 Fligner-Killeen test of homogeneity of variances
#> SE_2_+_32_35_37_38_ALDOA        Fligner-Killeen test of homogeneity of variances
#> SE_3_+_32_35_37_38_ANXA6        Fligner-Killeen test of homogeneity of variances
#> MXE_2_+_32_35_37_38_40_42_ALDOA Fligner-Killeen test of homogeneity of variances
#> MXE_3_+_32_35_37_38_40_42_ANXA6 Fligner-Killeen test of homogeneity of variances
#>                                 Fligner-Killeen data name Variance (Normal)
#> SE_1_+_32_35_37_38_ACTN1                 vector and group       0.000000000
#> MXE_1_+_32_35_37_38_40_42_ACTN1          vector and group       0.000000000
#> SE_2_+_32_35_37_38_ALDOA                 vector and group       0.026663118
#> SE_3_+_32_35_37_38_ANXA6                 vector and group       0.005240524
#> MXE_2_+_32_35_37_38_40_42_ALDOA          vector and group       0.019685695
#> MXE_3_+_32_35_37_38_40_42_ANXA6          vector and group       0.022583581
#>                                 Variance (Tumour) Median (Normal)
#> SE_1_+_32_35_37_38_ACTN1             0.0000000000       0.5000000
#> MXE_1_+_32_35_37_38_40_42_ACTN1      0.0000000000       0.5000000
#> SE_2_+_32_35_37_38_ALDOA             0.0027119499       0.7402597
#> SE_3_+_32_35_37_38_ANXA6             0.0003158887       0.4255319
#> MXE_2_+_32_35_37_38_40_42_ALDOA      0.0035641399       0.7500000
#> MXE_3_+_32_35_37_38_40_42_ANXA6      0.0027577471       0.2597403
#>                                 Median (Tumour) T-test statistic
#> SE_1_+_32_35_37_38_ACTN1              0.5000000               NA
#> MXE_1_+_32_35_37_38_40_42_ACTN1       0.5000000               NA
#> SE_2_+_32_35_37_38_ALDOA              0.6037736       1.04878918
#> SE_3_+_32_35_37_38_ANXA6              0.4418605      -0.02602006
#> MXE_2_+_32_35_37_38_40_42_ALDOA       0.5918367       1.73393229
#> MXE_3_+_32_35_37_38_40_42_ANXA6       0.3846154      -1.45806393
#>                                 T-test parameter T-test p-value
#> SE_1_+_32_35_37_38_ACTN1                      NA             NA
#> MXE_1_+_32_35_37_38_40_42_ACTN1               NA             NA
#> SE_2_+_32_35_37_38_ALDOA                2.402681      0.3881661
#> SE_3_+_32_35_37_38_ANXA6                2.240239      0.9813765
#> MXE_2_+_32_35_37_38_40_42_ALDOA         2.701223      0.1913725
#> MXE_3_+_32_35_37_38_40_42_ANXA6         2.481275      0.2586925
#>                                 T-test p-value (BH adjusted) T-test conf int1
#> SE_1_+_32_35_37_38_ACTN1                                  NA               NA
#> MXE_1_+_32_35_37_38_40_42_ACTN1                           NA               NA
#> SE_2_+_32_35_37_38_ALDOA                           0.5175548       -0.2604062
#> SE_3_+_32_35_37_38_ANXA6                           0.9813765       -0.1685048
#> MXE_2_+_32_35_37_38_40_42_ALDOA                    0.5173850       -0.1458915
#> MXE_3_+_32_35_37_38_40_42_ANXA6                    0.5173850       -0.4643286
#>                                 T-test conf int2 T-test estimate1
#> SE_1_+_32_35_37_38_ACTN1                      NA               NA
#> MXE_1_+_32_35_37_38_40_42_ACTN1               NA               NA
#> SE_2_+_32_35_37_38_ALDOA               0.4679678        0.7217532
#> SE_3_+_32_35_37_38_ANXA6               0.1662652        0.4353827
#> MXE_2_+_32_35_37_38_40_42_ALDOA        0.4511807        0.7434609
#> MXE_3_+_32_35_37_38_40_42_ANXA6        0.1963126        0.2671886
#>                                 T-test estimate2 T-test null value
#> SE_1_+_32_35_37_38_ACTN1                      NA                NA
#> MXE_1_+_32_35_37_38_40_42_ACTN1               NA                NA
#> SE_2_+_32_35_37_38_ALDOA               0.6179724                 0
#> SE_3_+_32_35_37_38_ANXA6               0.4365025                 0
#> MXE_2_+_32_35_37_38_40_42_ALDOA        0.5908163                 0
#> MXE_3_+_32_35_37_38_40_42_ANXA6        0.4011966                 0
#>                                 T-test stderr T-test alternative
#> SE_1_+_32_35_37_38_ACTN1                   NA               <NA>
#> MXE_1_+_32_35_37_38_40_42_ACTN1            NA               <NA>
#> SE_2_+_32_35_37_38_ALDOA           0.09895296          two.sided
#> SE_3_+_32_35_37_38_ANXA6           0.04303647          two.sided
#> MXE_2_+_32_35_37_38_40_42_ALDOA    0.08803377          two.sided
#> MXE_3_+_32_35_37_38_40_42_ANXA6    0.09190816          two.sided
#>                                           T-test method
#> SE_1_+_32_35_37_38_ACTN1                           <NA>
#> MXE_1_+_32_35_37_38_40_42_ACTN1                    <NA>
#> SE_2_+_32_35_37_38_ALDOA        Welch Two Sample t-test
#> SE_3_+_32_35_37_38_ANXA6        Welch Two Sample t-test
#> MXE_2_+_32_35_37_38_40_42_ALDOA Welch Two Sample t-test
#> MXE_3_+_32_35_37_38_40_42_ANXA6 Welch Two Sample t-test
#>                                                     T-test data name
#> SE_1_+_32_35_37_38_ACTN1                                        <NA>
#> MXE_1_+_32_35_37_38_40_42_ACTN1                                 <NA>
#> SE_2_+_32_35_37_38_ALDOA        vector[typeOne] and vector[!typeOne]
#> SE_3_+_32_35_37_38_ANXA6        vector[typeOne] and vector[!typeOne]
#> MXE_2_+_32_35_37_38_40_42_ALDOA vector[typeOne] and vector[!typeOne]
#> MXE_3_+_32_35_37_38_40_42_ANXA6 vector[typeOne] and vector[!typeOne]
#>                                  ∆ Variance    ∆ Median
#> SE_1_+_32_35_37_38_ACTN1        0.000000000  0.00000000
#> MXE_1_+_32_35_37_38_40_42_ACTN1 0.000000000  0.00000000
#> SE_2_+_32_35_37_38_ALDOA        0.023951168  0.13648616
#> SE_3_+_32_35_37_38_ANXA6        0.004924635 -0.01632855
#> MXE_2_+_32_35_37_38_40_42_ALDOA 0.016121555  0.15816327
#> MXE_3_+_32_35_37_38_40_42_ANXA6 0.019825834 -0.12487512
```
