# Multiple independence tests between reference groups and list of groups

Test multiple contingency tables comprised by two groups (one reference
group and another containing remaining elements) and provided groups.

## Usage

``` r
testGroupIndependence(ref, groups, elements, pvalueAdjust = "BH")
```

## Arguments

- ref:

  List of character: list of groups where each element contains the
  identifiers of respective elements

- groups:

  List of characters: list of groups where each element contains the
  identifiers of respective elements

- elements:

  Character: all available elements (if a data frame is given, its
  rownames will be used)

- pvalueAdjust:

  Character: method used to adjust p-values (see Details)

## Value

`multiGroupIndependenceTest` object, a data frame containing:

- attribute:

  Name of the original groups compared against the reference groups

- table:

  Contingency table used for testing

- pvalue:

  Fisher's exact test's p-value

## Details

The following methods for p-value adjustment are supported by using the
respective string in the `pvalueAdjust` argument:

- `none`: Do not adjust p-values

- `BH`: Benjamini-Hochberg's method (false discovery rate)

- `BY`: Benjamini-Yekutieli's method (false discovery rate)

- `bonferroni`: Bonferroni correction (family-wise error rate)

- `holm`: Holm's method (family-wise error rate)

- `hochberg`: Hochberg's method (family-wise error rate)

- `hommel`: Hommel's method (family-wise error rate)

## See also

[`parseCategoricalGroups()`](https://nuno-agostinho.github.io/psichomics/reference/parseCategoricalGroups.md)
and
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md)

Other functions for data grouping:
[`createGroupByAttribute()`](https://nuno-agostinho.github.io/psichomics/reference/createGroupByAttribute.md),
[`getGeneList()`](https://nuno-agostinho.github.io/psichomics/reference/getGeneList.md),
[`getSampleFromSubject()`](https://nuno-agostinho.github.io/psichomics/reference/getSampleFromSubject.md),
[`getSubjectFromSample()`](https://nuno-agostinho.github.io/psichomics/reference/getSubjectFromSample.md),
[`groupPerElem()`](https://nuno-agostinho.github.io/psichomics/reference/groupPerElem.md),
[`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md)

## Examples

``` r
elements <- paste("subjects", 1:10)
ref      <- elements[5:10]
groups   <- list(race=list(asian=elements[1:3],
                           white=elements[4:7],
                           black=elements[8:10]),
                 region=list(european=elements[c(4, 5, 9)],
                             african=elements[c(6:8, 10)]))
groupTesting <- testGroupIndependence(ref, groups, elements)
# View(groupTesting)
```
