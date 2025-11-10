# Process survival data to calculate survival curves

Process survival data to calculate survival curves

## Usage

``` r
processSurvData(
  event,
  timeStart,
  timeStop,
  followup,
  group,
  clinical,
  survTime = NULL
)
```

## Arguments

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

- group:

  Character: group relative to each subject

- clinical:

  Data frame: clinical data

- survTime:

  `survTime` object: Times to follow up, time start, time stop and event
  (optional)

## Value

Data frame with terms needed to calculate survival curves

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
