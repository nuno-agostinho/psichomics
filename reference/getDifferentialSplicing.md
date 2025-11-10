# Get or set differential splicing' elements for a data category

Get or set differential splicing' elements for a data category

## Usage

``` r
getDifferentialSplicing(category = getCategory())

setDifferentialSplicing(differential, category = getCategory())

getDifferentialSplicingFiltered(category = getCategory())

setDifferentialSplicingFiltered(differential, category = getCategory())

getDifferentialSplicingSurvival(category = getCategory())

setDifferentialSplicingSurvival(survival, category = getCategory())

getDifferentialSplicingResetPaging(category = getCategory())

setDifferentialSplicingResetPaging(reset, category = getCategory())

getDifferentialSplicingColumns(category = getCategory())

setDifferentialSplicingColumns(columns, category = getCategory())
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
[`getDifferentialExpression()`](https://nuno-agostinho.github.io/psichomics/reference/getDifferentialExpression.md),
[`getGlobal()`](https://nuno-agostinho.github.io/psichomics/reference/getGlobal.md),
[`getGroups()`](https://nuno-agostinho.github.io/psichomics/reference/getGroups.md),
[`getHighlightedPoints()`](https://nuno-agostinho.github.io/psichomics/reference/getHighlightedPoints.md),
[`getSelectedDataPanel()`](https://nuno-agostinho.github.io/psichomics/reference/getSelectedDataPanel.md)
