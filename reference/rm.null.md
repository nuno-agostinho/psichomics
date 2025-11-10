# Filter `NULL` elements from a vector or a list

Filter `NULL` elements from a vector or a list

## Usage

``` r
rm.null(v)
```

## Arguments

- v:

  Vector or list

## Value

Filtered vector or list with no `NULL` elements; if `v` is a vector
composed of `NULL` elements, returns a `NULL`; if `v` is a list of
`NULL` elements, returns an empty list
