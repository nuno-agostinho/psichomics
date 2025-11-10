# Get or set differential expression' elements for a data category

Get or set differential expression' elements for a data category

## Usage

``` r
getDifferentialExpression(category = getCategory())

setDifferentialExpression(differential, category = getCategory())

getDifferentialExpressionFiltered(category = getCategory())

setDifferentialExpressionFiltered(differential, category = getCategory())

getDifferentialExpressionSurvival(category = getCategory())

setDifferentialExpressionSurvival(survival, category = getCategory())

getDifferentialExpressionResetPaging(category = getCategory())

setDifferentialExpressionResetPaging(reset, category = getCategory())

getDifferentialExpressionColumns(category = getCategory())

setDifferentialExpressionColumns(columns, category = getCategory())
```

## Arguments

- category:

  Character: data category

- differential:

  Data frame or matrix: differential analyses table

- survival:

  Data frame or matrix: differential analyses' survival data

- reset:

  Character: reset paging of differential analyses table?

- columns:

  Character: differential analyses' column names

## Value

Getters return globally accessible data, whereas setters return `NULL`
as they are only used to modify the Shiny session's state

## Note

Needs to be called inside a reactive function

## See also

Other functions to get and set global variables:
[`getClinicalMatchFrom()`](https://nuno-agostinho.github.io/psichomics/reference/getClinicalMatchFrom.md),
[`getDifferentialSplicing()`](https://nuno-agostinho.github.io/psichomics/reference/getDifferentialSplicing.md),
[`getGlobal()`](https://nuno-agostinho.github.io/psichomics/reference/getGlobal.md),
[`getGroups()`](https://nuno-agostinho.github.io/psichomics/reference/getGroups.md),
[`getHighlightedPoints()`](https://nuno-agostinho.github.io/psichomics/reference/getHighlightedPoints.md),
[`getSelectedDataPanel()`](https://nuno-agostinho.github.io/psichomics/reference/getSelectedDataPanel.md)
