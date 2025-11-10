# Logic set to create an interactive [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

Logic set to create an interactive
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Usage

``` r
ggplotServer(
  input,
  output,
  id,
  plot = NULL,
  df = NULL,
  x = NULL,
  y = NULL,
  eventData = NULL
)

ggplotAuxServer(input, output, id)
```

## Arguments

- input:

  Shiny input

- output:

  Shiny output

- id:

  Character: identifier

- plot:

  Character: plot expression (if `NULL`, no plots are rendered)

- df:

  Data frame

- x:

  Character: name of the variable used for the X axis

- y:

  Character: name of the variable used for the Y axis

- eventData:

  Alternative splicing event information (if available)

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## Note

Insert `ggplotAuxSet` outside any observer (so it is only run once)
