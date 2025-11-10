# Parse categorical columns in a data frame

Retrieve elements grouped by their unique group based on each
categorical column

## Usage

``` r
parseCategoricalGroups(df)
```

## Arguments

- df:

  Data frame

## Value

List of lists containing values based on rownames of `df`

## See also

[`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)
and
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md)

## Examples

``` r
df <- data.frame("race"=c("caucasian", "caucasian", "asian"),
                 "gender"=c("male", "female", "male"))
rownames(df) <- paste("subject", 1:3)
parseCategoricalGroups(df)
#> named list()
```
