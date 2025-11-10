# Calculate optimal data cutoff that best separates survival curves

Uses [`stats::optim`](https://rdrr.io/r/stats/optim.html) with the Brent
method to test multiple cutoffs and to find the minimum log-rank
p-value.

## Usage

``` r
optimalSurvivalCutoff(
  clinical,
  data,
  censoring,
  event,
  timeStart,
  timeStop = NULL,
  followup = "days_to_last_followup",
  session = NULL,
  filter = TRUE,
  survTime = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- clinical:

  Data frame: clinical data

- data:

  Numeric: data values

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

- followup:

  Character: name of column containing follow up time

- session:

  Shiny session (only used for the visual interface)

- filter:

  Boolean or numeric: elements to use (all are used by default)

- survTime:

  `survTime` object: times to follow up, time start, time stop and event
  (optional)

- lower, upper:

  Bounds in which to search (if `NULL`, bounds are set to `lower = 0`
  and `upper = 1` if all data values are within that interval;
  otherwise, `lower = min(data, na.rm = TRUE)` and
  `upper = max(data, na.rm = TRUE)`)

## Value

List containing the optimal cutoff (`par`) and the corresponding p-value
(`value`)

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
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

psi <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
opt <- optimalSurvivalCutoff(clinical, psi, "right", event, timeStart)
```
