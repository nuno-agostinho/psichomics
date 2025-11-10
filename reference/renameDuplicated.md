# Rename vector to avoid duplicated values with another vector

Renames values by adding an index to the end of duplicates. This allows
to prepare unique values in two vectors before a merge, for instance.

## Usage

``` r
renameDuplicated(check, comp)
```

## Arguments

- check:

  Character: values to rename if duplicated

- comp:

  Character: values to compare with

## Value

Character vector with renamed values if duplicated; else, it returns the
usual values. It does not return the comparator values.

## Examples

``` r
psichomics:::renameDuplicated(check = c("blue", "red"), comp = c("green",
                                                                 "blue"))
#> [1] "blue (1)" "red"     
```
