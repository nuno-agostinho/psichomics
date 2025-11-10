# Create density sparklines for inclusion levels

Create density sparklines for inclusion levels

## Usage

``` r
createDensitySparklines(
  data,
  events,
  areSplicingEvents = TRUE,
  groups = NULL,
  geneExpr = NULL,
  inputID = "sparklineInput"
)
```

## Arguments

- data:

  Character: HTML-formatted data series of interest

- events:

  Character: event identifiers

- areSplicingEvents:

  Boolean: are these splicing events (TRUE) or gene expression (FALSE)?

- groups:

  Character: name of the groups used for differential analyses

- geneExpr:

  Character: name of the gene expression dataset

- inputID:

  Character: identifier of input to get attributes of clicked event
  (Shiny only)

## Value

HTML element with sparkline data
