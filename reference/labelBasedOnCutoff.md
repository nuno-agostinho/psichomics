# Label groups based on a given cutoff

Label groups based on a given cutoff

## Usage

``` r
labelBasedOnCutoff(data, cutoff, label = NULL, gte = TRUE)
```

## Arguments

- data:

  Numeric: test data

- cutoff:

  Numeric: test cutoff

- label:

  Character: label to prefix group names

- gte:

  Boolean: test using greater than or equal than cutoff (`TRUE`) or less
  than or equal than cutoff (`FALSE`)?

## Value

Labelled groups

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
[`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md),
[`survfit.survTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md),
[`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)

## Examples

``` r
labelBasedOnCutoff(data=c(1, 0, 0, 1, 0, 1), cutoff=0.5)
#> [1] "&gt;= 0.5" "&lt; 0.5"  "&lt; 0.5"  "&gt;= 0.5" "&lt; 0.5"  "&gt;= 0.5"

labelBasedOnCutoff(data=c(1, 0, 0, 1, 0, 1), cutoff=0.5, "Ratio")
#> [1] "Ratio &gt;= 0.5" "Ratio &lt; 0.5"  "Ratio &lt; 0.5"  "Ratio &gt;= 0.5"
#> [5] "Ratio &lt; 0.5"  "Ratio &gt;= 0.5"

# Use "greater than" instead of "greater than or equal to"
labelBasedOnCutoff(data=c(1, 0, 0, 0.5, 0, 1), cutoff=0.5, gte=FALSE)
#> [1] "&gt; 0.5"  "&lt;= 0.5" "&lt;= 0.5" "&lt;= 0.5" "&lt;= 0.5" "&gt; 0.5" 
```
