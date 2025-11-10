# Get or set groups

Get or set groups

## Usage

``` r
getGroups(
  type = c("Patients", "Samples", "ASevents", "Genes"),
  complete = FALSE,
  category = getCategory()
)

setGroups(
  type = c("Patients", "Samples", "ASevents", "Genes"),
  groups,
  category = getCategory()
)
```

## Arguments

- type:

  Character: type of groups (either `Patients`, `Samples`, `ASevents` or
  `Genes`)

- complete:

  Boolean: return all the information on groups (`TRUE`) or just the
  group names and respective indexes (`FALSE`)?

- category:

  Character: data category

- groups:

  Matrix: groups of dataset

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
[`getHighlightedPoints()`](https://nuno-agostinho.github.io/psichomics/reference/getHighlightedPoints.md),
[`getSelectedDataPanel()`](https://nuno-agostinho.github.io/psichomics/reference/getSelectedDataPanel.md)
