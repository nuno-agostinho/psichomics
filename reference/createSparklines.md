# Create sparkline charts to be used in a data table

Create sparkline charts to be used in a data table

## Usage

``` r
createSparklines(
  hc,
  data,
  events,
  groups = NULL,
  geneExpr = NULL,
  inputID = "sparklineInput",
  ...
)
```

## Arguments

- hc:

  `highchart` object

- data:

  Character: HTML-formatted data series of interest

- events:

  Character: event identifiers

- groups:

  Character: name of the groups used for differential analyses

- geneExpr:

  Character: name of the gene expression dataset

- inputID:

  Character: identifier of input to get attributes of clicked event
  (Shiny only)

- id:

  Character: Shiny input identifier

## Value

HTML element with sparkline data
