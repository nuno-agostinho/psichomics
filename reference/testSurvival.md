# Test the survival difference between groups of subjects

Test the survival difference between groups of subjects

## Usage

``` r
testSurvival(survTerms, ...)
```

## Arguments

- survTerms:

  `survTerms` object: survival terms obtained after running
  `processSurvTerms` (see examples)

- ...:

  Arguments passed on to
  [`survival::survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html)

  `subset`

  :   expression indicating which subset of the rows of data should be
      used in the fit. This can be a logical vector (which is replicated
      to have length equal to the number of observations), a numeric
      vector indicating which observation numbers are to be included (or
      excluded if negative), or a character vector of row names to be
      included. All observations are included by default.

  `na.action`

  :   a missing-data filter function. This is applied to the
      `model.frame` after any subset argument has been used. Default is
      `options()$na.action`.

  `rho`

  :   a scalar parameter that controls the type of test.

  `timefix`

  :   process times through the `aeqSurv` function to eliminate
      potential roundoff issues.

## Value

p-value of the survival difference or `NA`

## Note

Instead of raising errors, returns `NA`

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
[`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md),
[`survfit.survTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md)

## Examples

``` r
require("survival")
data <- aml
timeStart  <- "event"
event      <- "event"
followup   <- "time"
data$event <- NA
data$event[aml$status == 1] <- aml$time[aml$status == 1]
censoring  <- "right"
formulaStr <- "x"
survTerms <- processSurvTerms(data, censoring=censoring, event=event,
                              timeStart=timeStart, followup=followup,
                              formulaStr=formulaStr)
testSurvival(survTerms)
#> [1] 0.0653
```
