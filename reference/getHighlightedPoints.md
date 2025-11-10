# Get or set points or regions for plots

Get or set points or regions for plots

## Usage

``` r
getHighlightedPoints(id, category = getCategory())

setHighlightedPoints(id, events, category = getCategory())

getZoom(id, category = getCategory())

setZoom(id, zoom, category = getCategory())

getSelectedPoints(id, category = getCategory())

setSelectedPoints(id, events, category = getCategory())

getLabelledPoints(id, category = getCategory())

setLabelledPoints(id, events, category = getCategory())
```

## Arguments

- id:

  Character: identifier

- category:

  Character: data category

- events:

  Integer: index of events

- zoom:

  Integer: range of X and Y coordinates for zooming

## Value

Getters return globally accessible data, whereas setters return `NULL`
as they are only used to modify the Shiny session's state

## Note

Needs to be called inside a reactive function

## See also

Other functions to get and set global variables:
[`getClinicalMatchFrom()`](https://nuno-agostinho.github.io/psichomics/reference/getClinicalMatchFrom.md),
[`getDifferentialExpression()`](https://nuno-agostinho.github.io/psichomics/reference/getDifferentialExpression.md),
[`getDifferentialSplicing()`](https://nuno-agostinho.github.io/psichomics/reference/getDifferentialSplicing.md),
[`getGlobal()`](https://nuno-agostinho.github.io/psichomics/reference/getGlobal.md),
[`getGroups()`](https://nuno-agostinho.github.io/psichomics/reference/getGroups.md),
[`getSelectedDataPanel()`](https://nuno-agostinho.github.io/psichomics/reference/getSelectedDataPanel.md)
