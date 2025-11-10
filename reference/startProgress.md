# Create, set and terminate a progress object

Create, set and terminate a progress object

## Usage

``` r
startProgress(
  message,
  divisions,
  global = if (isRunning()) sharedData else getHidden()
)

updateProgress(
  message = "Loading...",
  value = NULL,
  max = NULL,
  detail = NULL,
  divisions = NULL,
  global = if (isRunning()) sharedData else getHidden(),
  console = TRUE
)

closeProgress(
  message = NULL,
  global = if (isRunning()) sharedData else getHidden()
)
```

## Arguments

- message:

  Character: progress message

- divisions:

  Integer: number of divisions in the progress bar

- global:

  Shiny's global variable

- value:

  Integer: current progress value

- max:

  Integer: maximum progress value

- detail:

  Character: detailed message

- console:

  Boolean: print message to console?

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## Details

If `divisions` is not `NULL`, a progress bar starts with the given
divisions. If `value = NULL`, the progress bar increments one unit;
otherwise, the progress bar increments `value`.
