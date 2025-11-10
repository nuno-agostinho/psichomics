# Process survival curves terms to calculate survival curves

Process survival curves terms to calculate survival curves

## Usage

``` r
processSurvTerms(
  clinical,
  censoring,
  event,
  timeStart,
  timeStop = NULL,
  group = NULL,
  formulaStr = NULL,
  coxph = FALSE,
  scale = "days",
  followup = "days_to_last_followup",
  survTime = NULL
)
```

## Arguments

- clinical:

  Data frame: clinical data

- censoring:

  Character: censor using `left`, `right`, `interval` or `interval2`

- event:

  Character: name of column containing time of the event of interest

- timeStart:

  Character: name of column containing starting time of the interval or
  follow up time

- timeStop:

  Character: name of column containing ending time of the interval (only
  relevant for interval censoring)

- group:

  Character: group relative to each subject

- formulaStr:

  Character: formula to use

- coxph:

  Boolean: fit a Cox proportional hazards regression model?

- scale:

  Character: rescale the survival time to `days`, `weeks`, `months` or
  `years`

- followup:

  Character: name of column containing follow up time

- survTime:

  `survTime` object: times to follow up, time start, time stop and event
  (optional)

## Value

A list with a `formula` object and a data frame with terms needed to
calculate survival curves

## Details

The `event` time is only used to determine whether the event has
occurred (`1`) or not (`0`) in case of missing values.

If `survTime = NULL`, survival times are obtained from the clinical
dataset according to the names given in `timeStart`, `timeStop`, `event`
and `followup`. This may become quite slow when used in a loop. If the
aforementioned variables are constant, consider running
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md)
outside the loop and using its output via the `survTime` argument of
this function (see Examples).

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md),
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

# If running multiple times, consider calculating survTime only once
survTime <- getAttributesTime(clinical, event, timeStart)
for (i in seq(5)) {
  survTerms <- processSurvTerms(clinical, censoring="right", event,
                                timeStart, formulaStr=formulaStr,
                                survTime=survTime)
}
```
