# Set the status of a process to style a given button

- `startProcess`: Style button to show a process is in progress

- `endProcess`: Style button to show a process finished; also, closes
  the progress bar (if `closeProgressbar = TRUE`) and prints the
  difference between the current time and `time`

## Usage

``` r
startProcess(id)

endProcess(id, time = NULL, closeProgressBar = TRUE)
```

## Arguments

- id:

  Character: button identifier

- time:

  `POSIXct` object: start time needed to show the interval time (if
  `NULL`, the time interval is not displayed)

- closeProgressBar:

  Boolean: close progress bar?

## Value

`startProcess` returns the start time of the process (may be used as the
`time` argument to `endProcess`), whereas `endProcess` returns the
difference between current time and `time` (or `NULL` if `time` is not
specified)
