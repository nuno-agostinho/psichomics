# Create the interface for the tooltip of a plot

Create the interface for the tooltip of a plot

## Usage

``` r
ggplotTooltip(df, hover, x, y, eventData = NULL)
```

## Arguments

- df:

  Data frame

- hover:

  Mouse hover information for a given plot as retrieved from
  [`hoverOpts`](https://rdrr.io/pkg/shiny/man/clickOpts.html)

- x:

  Character: name of the variable used for the X axis

- y:

  Character: name of the variable used for the Y axis

- eventData:

  Alternative splicing event information (if available)

## Value

HTML elements
