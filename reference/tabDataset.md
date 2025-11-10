# Creates a `tabPanel` template for a `datatable` with a title and description

Creates a `tabPanel` template for a `datatable` with a title and
description

## Usage

``` r
tabDataset(
  ns,
  title,
  tableId,
  columns,
  visCols,
  data,
  description = NULL,
  icon = NULL
)
```

## Arguments

- ns:

  Namespace function

- title:

  Character: tab title

- tableId:

  Character: id of the `datatable`

- columns:

  Character: column names of the `datatable`

- visCols:

  Boolean: visible columns

- data:

  Data frame: dataset of interest

- description:

  Character: description of the table (optional)

- icon:

  Character: list containing an item named `symbol` (FontAwesome icon
  name) and another one named `colour` (background colour)

## Value

HTML elements
