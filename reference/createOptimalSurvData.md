# Create survival data based on a PSI cutoff

Data is presented in the table for statistical analyses

## Usage

``` r
createOptimalSurvData(
  eventPSI,
  clinical,
  censoring,
  event,
  timeStart,
  timeStop,
  match,
  patients,
  samples
)
```

## Arguments

- eventPSI:

  Numeric: alternative splicing quantification for multiple samples
  relative to a single splicing event

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

- match:

  Matrix: match between samples and subjects

- patients:

  Character: subject identifiers (only required if the `clinical`
  argument is not handed)

- samples:

  Character: samples to use when assigning values per subject (if
  `NULL`, all samples will be used)

## Value

Survival data including optimal PSI cutoff, minimal survival p-value and
HTML element required to plot survival curves
