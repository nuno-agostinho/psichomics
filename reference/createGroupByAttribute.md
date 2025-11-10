# Split elements into groups based on a given column of a dataset

Elements are identified by their respective row name.

## Usage

``` r
createGroupByAttribute(col, dataset)
```

## Arguments

- col:

  Character: column name

- dataset:

  Matrix or data frame: dataset

## Value

Named list with each unique value from a given column and respective
elements

## See also

Other functions for data grouping:
[`getGeneList()`](https://nuno-agostinho.github.io/psichomics/reference/getGeneList.md),
[`getSampleFromSubject()`](https://nuno-agostinho.github.io/psichomics/reference/getSampleFromSubject.md),
[`getSubjectFromSample()`](https://nuno-agostinho.github.io/psichomics/reference/getSubjectFromSample.md),
[`groupPerElem()`](https://nuno-agostinho.github.io/psichomics/reference/groupPerElem.md),
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md),
[`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)

## Examples

``` r
df <- data.frame(gender=c("male", "female"),
                 stage=paste("stage", c(1, 3, 1, 4, 2, 3, 2, 2)))
rownames(df) <- paste0("subject-", LETTERS[1:8])
createGroupByAttribute(col="stage", dataset=df)
#> $`stage 1`
#> [1] "subject-A" "subject-C"
#> 
#> $`stage 2`
#> [1] "subject-E" "subject-G" "subject-H"
#> 
#> $`stage 3`
#> [1] "subject-B" "subject-F"
#> 
#> $`stage 4`
#> [1] "subject-D"
#> 
```
