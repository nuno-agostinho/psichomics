# Prepare to create group according to specific details

Prepare to create group according to specific details

## Usage

``` r
createGroup(
  session,
  input,
  output,
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

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
