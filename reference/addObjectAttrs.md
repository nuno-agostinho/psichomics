# Set attributes to an object

Set attributes to an object

## Usage

``` r
addObjectAttrs(object, ..., replace = TRUE)
```

## Arguments

- object:

  Object

- ...:

  Named parameters to convert to attributes

- replace:

  Boolean: replace an attribute if already set?

## Value

Object with attributes set

## Examples

``` r
ll <- list(a="hey", b="there")
psichomics:::addObjectAttrs(ll, "words"=2, "language"="English")
#> $a
#> [1] "hey"
#> 
#> $b
#> [1] "there"
#> 
#> attr(,"words")
#> [1] 2
#> attr(,"language")
#> [1] "English"
```
