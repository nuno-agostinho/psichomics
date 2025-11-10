# Assign one group to each element

Assign one group to each element

## Usage

``` r
groupPerElem(groups, elem = NULL, outerGroupName = NA)
```

## Arguments

- groups:

  List of integers: groups of elements

- elem:

  Character: all elements available

- outerGroupName:

  Character: name to give to outer group (if `NULL`, only show elements
  matched to their respective groups)

## Value

Character vector where each element corresponds to the group of the
respective element

## See also

Other functions for data grouping:
[`createGroupByAttribute()`](https://nuno-agostinho.github.io/psichomics/reference/createGroupByAttribute.md),
[`getGeneList()`](https://nuno-agostinho.github.io/psichomics/reference/getGeneList.md),
[`getSampleFromSubject()`](https://nuno-agostinho.github.io/psichomics/reference/getSampleFromSubject.md),
[`getSubjectFromSample()`](https://nuno-agostinho.github.io/psichomics/reference/getSubjectFromSample.md),
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md),
[`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)

## Examples

``` r
groups <- list(1:3, 4:7, 8:10)
names(groups) <- paste("Stage", 1:3)
groupPerElem(groups)
#>         1         2         3         4         5         6         7         8 
#> "Stage 1" "Stage 1" "Stage 1" "Stage 2" "Stage 2" "Stage 2" "Stage 2" "Stage 3" 
#>         9        10 
#> "Stage 3" "Stage 3" 
```
