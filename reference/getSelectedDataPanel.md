# Get or set selected panel in data section

Get or set selected panel in data section

## Usage

``` r
getSelectedDataPanel()

setSelectedDataPanel(id)
```

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
[`getHighlightedPoints()`](https://nuno-agostinho.github.io/psichomics/reference/getHighlightedPoints.md)
