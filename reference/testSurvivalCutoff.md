# Test the survival difference between two survival groups given a cutoff

Test the survival difference between two survival groups given a cutoff

## Usage

``` r
testSurvivalCutoff(
  cutoff,
  data,
  filter = TRUE,
  clinical,
  ...,
  session = NULL,
  survivalInfo = FALSE
)
```

## Arguments

- cutoff:

  Numeric: Cutoff of interest

- data:

  Numeric: elements of interest to test against the cutoff

- filter:

  Boolean or numeric: elements to use (all are used by default)

- clinical:

  Data frame: clinical data

- ...:

  Arguments passed on to
  [`processSurvTerms`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md)

  `censoring`

  :   Character: censor using `left`, `right`, `interval` or `interval2`

  `scale`

  :   Character: rescale the survival time to `days`, `weeks`, `months`
      or `years`

  `formulaStr`

  :   Character: formula to use

  `coxph`

  :   Boolean: fit a Cox proportional hazards regression model?

  `survTime`

  :   `survTime` object: times to follow up, time start, time stop and
      event (optional)

  `event`

  :   Character: name of column containing time of the event of interest

  `timeStart`

  :   Character: name of column containing starting time of the interval
      or follow up time

  `timeStop`

  :   Character: name of column containing ending time of the interval
      (only relevant for interval censoring)

  `followup`

  :   Character: name of column containing follow up time

- session:

  Shiny session

- survivalInfo:

  Boolean: return extra survival information

## Value

p-value of the survival difference
