# Create scatter plot

Create a scatter plot using `highcharter`

## Usage

``` r
hc_scatter(
  hc,
  x,
  y,
  z = NULL,
  label = NULL,
  showInLegend = FALSE,
  color = NULL,
  ...
)
```

## Arguments

- hc:

  `Highchart` object

- x:

  Numeric: X axis

- y:

  Numeric: Y axis

- z:

  Numeric: Z axis to set the bubble size (optional)

- label:

  Character: data label for each point (optional)

- showInLegend:

  Boolean: show the data in the legend box?

- color:

  Character: series colour

- ...:

  Arguments passed on to
  [`highcharter::hc_add_series`](https://jkunst.com/highcharter/reference/hc_add_series.html)

  :   

## Value

`highcharter` object containing information for a scatter plot
