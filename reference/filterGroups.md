# Filter groups with less data points than the threshold

Groups containing a number of non-missing values less than the threshold
are discarded.

## Usage

``` r
filterGroups(vector, group, threshold = 1)
```

## Arguments

- vector:

  Character: elements

- group:

  Character: respective group of each elements

- threshold:

  Integer: number of valid non-missing values by group

## Value

Named vector with filtered elements from valid groups. The group of the
respective element is given as an attribute.

## Examples

``` r
# Removes groups with less than two elements
vec <- 1:6
names(vec) <- paste("sample", letters[1:6])
filterGroups(vec, c("A", "B", "B", "C", "D", "D"), threshold=2)
#> sample b sample c sample e sample f 
#>        2        3        5        6 
#> attr(,"Groups")
#> [1] "B" "B" "D" "D"
```
