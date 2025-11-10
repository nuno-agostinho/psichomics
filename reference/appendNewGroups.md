# Append new groups to already existing groups

Retrieve previous groups, rename duplicated group names in the new
groups and append new groups to the previous ones

## Usage

``` r
appendNewGroups(type, new, clearOld = FALSE)
```

## Arguments

- type:

  Character: type of groups (either `Patients`, `Samples`, `ASevents` or
  `Genes`)

- new:

  Rows of groups to be added

- clearOld:

  Boolean: clear old groups?

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
