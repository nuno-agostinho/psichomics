# Set new groups according to the user input

Set new groups according to the user input

## Usage

``` r
createGroupFromInput(
  session,
  input,
  output,
  dataset,
  id,
  type,
  selected = NULL,
  expr = NULL,
  groupNames = NULL
)
```

## Arguments

- session:

  Shiny session

- input:

  Shiny input

- output:

  Shiny output

- dataset:

  Data frame or matrix: dataset of interest

- id:

  Character: identifier of the group selection

- type:

  Character: type of group to create

- selected:

  Character: selected item

- expr:

  Character: expression

- groupNames:

  Character: group names

## Value

Matrix with the group names and respective elements
