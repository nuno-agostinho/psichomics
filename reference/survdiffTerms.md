# Test Survival Curve Differences

Tests if there is a difference between two or more survival curves using
the \\G^\rho\\ family of tests, or for a single curve against a known
alternative.

## Usage

``` r
survdiffTerms(survTerms, ...)
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

`survfit` object. See `survfit.object` for details. Methods defined for
`survfit` objects are `print`, `plot`, `lines`, and `points`.

## Description

This function implements the G-rho family of Harrington and Fleming
(1982), with weights on each death of \\S(t)^\rho\\, where \\S(t)\\ is
the Kaplan-Meier estimate of survival. With `rho = 0` this is the
log-rank or Mantel-Haenszel test, and with `rho = 1` it is equivalent to
the Peto & Peto modification of the Gehan-Wilcoxon test.

Peto and Peto show that the Gehan-Wilcoxon test can be badly biased if
the two groups have different censoring patterns, and proposed an
alternative. Prentice and Marek later showed an actual example where
this issue occurs. For most data sets the Gehan-Wilcoxon and
Peto-Peto-Prentice variant will hardly differ, however.

If the right hand side of the formula consists only of an offset term,
then a one sample test is done. To cause missing values in the
predictors to be treated as a separate group, rather than being omitted,
use the `factor` function with its `exclude` argument to recode the
righ-hand-side covariate.

Note that the ordinary log-rank test is equivalent to the score test
from a Cox model, using the Breslow approximation for ties. Use the Cox
model form for more complex models, e.g., time-dependent covariates.

## References

Harrington, D. P. and Fleming, T. R. (1982). A class of rank test
procedures for censored survival data. Biometrika, 553-566.

Peto R. Peto and Peto, J. (1972) Asymptotically efficient rank invariant
test procedures (with discussion), JRSSA, 185-206.

Prentice, R. and Marek, P. (1979) A qualitative discrepancy between
censored data rank tests, Biometics, 861–867.

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
[`survfit.survTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md),
[`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)

## Examples

``` r
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
survdiffTerms(survTerms)
#> Call:
#> survdiff(formula = survTerms$form, data = survTerms$survTime)
#> 
#>                                                                 N Observed
#> patient.stage_event.pathologic_stage=i, patient.gender=female   1        0
#> patient.stage_event.pathologic_stage=ii, patient.gender=female  1        0
#> patient.stage_event.pathologic_stage=ii, patient.gender=male    1        1
#> patient.stage_event.pathologic_stage=iii, patient.gender=male   1        0
#> patient.stage_event.pathologic_stage=iv, patient.gender=female  1        1
#> patient.stage_event.pathologic_stage=iv, patient.gender=male    1        1
#>                                                                 Expected
#> patient.stage_event.pathologic_stage=i, patient.gender=female      0.167
#> patient.stage_event.pathologic_stage=ii, patient.gender=female     0.917
#> patient.stage_event.pathologic_stage=ii, patient.gender=male       0.917
#> patient.stage_event.pathologic_stage=iii, patient.gender=male      0.417
#> patient.stage_event.pathologic_stage=iv, patient.gender=female     0.167
#> patient.stage_event.pathologic_stage=iv, patient.gender=male       0.417
#>                                                                 (O-E)^2/E
#> patient.stage_event.pathologic_stage=i, patient.gender=female     0.16667
#> patient.stage_event.pathologic_stage=ii, patient.gender=female    0.91667
#> patient.stage_event.pathologic_stage=ii, patient.gender=male      0.00758
#> patient.stage_event.pathologic_stage=iii, patient.gender=male     0.41667
#> patient.stage_event.pathologic_stage=iv, patient.gender=female    4.16667
#> patient.stage_event.pathologic_stage=iv, patient.gender=male      0.81667
#>                                                                 (O-E)^2/V
#> patient.stage_event.pathologic_stage=i, patient.gender=female       0.200
#> patient.stage_event.pathologic_stage=ii, patient.gender=female      1.458
#> patient.stage_event.pathologic_stage=ii, patient.gender=male        0.012
#> patient.stage_event.pathologic_stage=iii, patient.gender=male       0.532
#> patient.stage_event.pathologic_stage=iv, patient.gender=female      5.000
#> patient.stage_event.pathologic_stage=iv, patient.gender=male        1.043
#> 
#>  Chisq= 7.3  on 5 degrees of freedom, p= 0.2 
```
