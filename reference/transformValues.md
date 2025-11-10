# Transform values as per a given type of transformation

Transform values as per a given type of transformation

## Usage

``` r
transformValues(val, type, avoidZero = TRUE)
```

## Arguments

- val:

  Integer: values to transform

- type:

  Character: type of transformation

- avoidZero:

  Boolean: add the smallest non-zero number available
  (`.Machine$double.xmin`) to avoid infinity values following
  log-transformation (may not be plotted); useful for p-values of 0

## Value

Integer containing transformed values
