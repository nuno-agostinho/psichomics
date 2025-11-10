# Create plot for events

Create plot for events

## Usage

``` r
createEventPlotting(
  df,
  x,
  y,
  params,
  highlightX,
  highlightY,
  highlightParams,
  selected,
  selectedParams,
  labelled,
  labelledParams,
  xlim,
  ylim
)
```

## Arguments

- df:

  Data frame

- x:

  Character: name of the variable used for the X axis

- y:

  Character: name of the variable used for the Y axis

- params:

  List of parameters to pass to
  [`geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  related to most points

- highlightX:

  Integer: region of points in X axis to highlight

- highlightY:

  Integer: region of points in Y axis to highlight

- highlightParams:

  List of parameters to pass to
  [`geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  related to highlighted points

- selected:

  Integer: index of rows/points to be coloured

- selectedParams:

  List of parameters to pass to
  [`geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  related to selected points

- labelled:

  Integer: index of rows/points to be labelled

- labelledParams:

  List of parameters to pass to
  [`ggrepel::geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
  related to labelled points

- xlim:

  Numeric: limits of X axis

- ylim:

  Numeric: limits of Y axis

## Value

List containing HTML elements and highlighted points
