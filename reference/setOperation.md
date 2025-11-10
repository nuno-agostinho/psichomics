# Perform set operations on selected groups

Perform set operations on selected groups

## Usage

``` r
setOperation(
  operation,
  groups,
  selected,
  symbol = " ",
  groupName = NULL,
  first = NULL,
  second = NULL,
  matches = NULL,
  type = "Samples",
  assignColoursToGroups = FALSE
)
```

## Arguments

- operation:

  Character: set operation

- groups:

  Matrix: groups

- selected:

  Integer: index of rows regarding selected groups

- symbol:

  Character: Unicode symbol to visually indicate the operation performed

- groupName:

  Character: group name (automatically created if `NULL` or `""`)

- first:

  Character: identifiers of the first element (required when performing
  the `complement` operation)

- second:

  Character: identifiers of the second element (required when performing
  the `complement` operation)

- matches:

  Character: match between samples (as names) and subjects (as values)

- type:

  Character: type of group where set operations are to be performed

- assignColoursToGroups:

  Boolean: assign colours to new groups?

## Value

Matrix containing groups (new group is in the first row)
