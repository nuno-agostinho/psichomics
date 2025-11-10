# Trims whitespace from a word

Trims whitespace from a word

## Usage

``` r
trimWhitespace(word)
```

## Arguments

- word:

  Character to trim

## Value

Character without whitespace

## Examples

``` r
psichomics:::trimWhitespace("    hey   there     ")
#> [1] "hey there"
psichomics:::trimWhitespace(c("pineapple    ", "one two three",
                              " sunken    ship   "))
#> [1] "pineapple"     "one two three" "sunken ship"  
```
