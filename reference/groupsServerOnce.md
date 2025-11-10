# Server function for data grouping (one call)

These functions only run once instead of running for every instance of
groups

## Usage

``` r
groupsServerOnce(input, output, session)
```

## Arguments

- input:

  Shiny input

- output:

  Shiny output

- session:

  Shiny session

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
