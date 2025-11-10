# Set operations on groups

This function can be used on groups to merge, intersect, subtract, etc.

## Usage

``` r
operateOnGroups(
  input,
  session,
  operation,
  buttonId,
  symbol = " ",
  type,
  sharedData = sharedData
)
```

## Arguments

- input:

  Shiny input

- session:

  Shiny session

- operation:

  Character: set operation

- buttonId:

  Character: ID of the button to trigger operation

- symbol:

  Character: Unicode symbol to visually indicate the operation performed

- type:

  Character: type of group where set operations are to be performed

- sharedData:

  Shiny app's global variable

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
