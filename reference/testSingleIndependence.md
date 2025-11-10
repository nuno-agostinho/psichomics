# Multiple independence tests between a reference group and list of groups

Uses Fisher's exact test.

## Usage

``` r
testSingleIndependence(ref, groups, elements, pvalueAdjust = "BH")
```

## Arguments

- ref:

  Character: identifier of elements in reference group

- groups:

  List of characters: list of groups where each element contains the
  identifiers of respective elements

- elements:

  Character: all subject identifiers

- pvalueAdjust:

  Character: method used to adjust p-values (see Details)

## Value

Returns a `groupIndependenceTest` object: a list where each element is a
list containing:

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
