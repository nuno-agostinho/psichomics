# Import groups from a file

Import groups from a file

## Usage

``` r
importGroupsFrom(
  file,
  uniqueElems = NULL,
  matchingElems = NULL,
  match = NULL,
  type = NULL
)
```

## Arguments

- file:

  Character: path to file

- uniqueElems:

  Character: vector of unique elements (samples or alternative splicing
  events)

- matchingElems:

  Character: vector of matching elements (subjects or genes)

- match:

  Match between elements within groups

## Value

Matrix with groups
