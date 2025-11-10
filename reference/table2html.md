# Create HTML table from data frame or matrix

Create HTML table from data frame or matrix

## Usage

``` r
table2html(
  data,
  rownames = TRUE,
  colnames = TRUE,
  class = NULL,
  style = NULL,
  thead = FALSE
)
```

## Arguments

- data:

  Data frame or matrix

- rownames:

  Boolean: print row names?

- colnames:

  Boolean: print column names?

- class:

  Character: table class

- style:

  Character: table style

- thead:

  Boolean: add a `thead` tag to the first row?

## Value

HTML elements
