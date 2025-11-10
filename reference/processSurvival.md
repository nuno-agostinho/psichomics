# Check if survival analyses successfully completed or returned errors

Check if survival analyses successfully completed or returned errors

## Usage

``` r
processSurvival(session, ...)
```

## Arguments

- session:

  Shiny session

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

  `group`

  :   Character: group relative to each subject

  `clinical`

  :   Data frame: clinical data

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

## Value

List with survival analysis results
