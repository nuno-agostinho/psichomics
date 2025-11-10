# Create survival curves

Create survival curves

## Usage

``` r
# S3 method for class 'survTerms'
survfit(formula, ...)
```

## Arguments

- formula:

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

`survfit` object. See `survfit.object` for details. Methods defined for
`survfit` objects are `print`, `plot`, `lines`, and `points`.

## Details

A survival curve is based on a tabulation of the number at risk and
number of events at each unique death time. When time is a floating
point number the definition of "unique" is subject to interpretation.
The code uses factor() to define the set. For further details see the
documentation for the appropriate method, i.e.,
[`?survfit.formula`](https://rdrr.io/pkg/survival/man/survfit.formula.html)
or
[`?survfit.coxph`](https://rdrr.io/pkg/survival/man/survfit.coxph.html).

A survfit object may contain a single curve, a set of curves (vector), a
matrix of curves, or even a 3 way array: `dim(fit)` will reveal the
dimensions. Predicted curves from a `coxph` model have one row for each
stratum in the Cox model fit and one column for each specified covariate
set. Curves from a multi-state model have one row for each stratum and a
column for each state, the strata correspond to predictors on the right
hand side of the equation. The default printing and plotting order for
curves is by column, as with other matrices.

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
[`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)

## Examples

``` r
library("survival")
clinical <- read.table(text = "2549   NA ii  female
                                840   NA i   female
                                 NA 1204 iv    male
                                 NA  383 iv  female
                               1293   NA iii   male
                                 NA 1355 ii    male")
names(clinical) <- c("patient.days_to_last_followup",
                     "patient.days_to_death",
                     "patient.stage_event.pathologic_stage",
                     "patient.gender")
timeStart  <- "days_to_death"
event      <- "days_to_death"
formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
                               formulaStr=formulaStr)
survfit(survTerms)
#> Call: survfit(formula = formula$form, data = formula$survTime)
#> 
#>                                                                 n events median
#> patient.stage_event.pathologic_stage=i, patient.gender=female   1      0     NA
#> patient.stage_event.pathologic_stage=ii, patient.gender=female  1      0     NA
#> patient.stage_event.pathologic_stage=ii, patient.gender=male    1      1   1355
#> patient.stage_event.pathologic_stage=iii, patient.gender=male   1      0     NA
#> patient.stage_event.pathologic_stage=iv, patient.gender=female  1      1    383
#> patient.stage_event.pathologic_stage=iv, patient.gender=male    1      1   1204
#>                                                                 0.95LCL 0.95UCL
#> patient.stage_event.pathologic_stage=i, patient.gender=female        NA      NA
#> patient.stage_event.pathologic_stage=ii, patient.gender=female       NA      NA
#> patient.stage_event.pathologic_stage=ii, patient.gender=male         NA      NA
#> patient.stage_event.pathologic_stage=iii, patient.gender=male        NA      NA
#> patient.stage_event.pathologic_stage=iv, patient.gender=female       NA      NA
#> patient.stage_event.pathologic_stage=iv, patient.gender=male         NA      NA
```
