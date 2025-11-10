# Blend two HEX colours

Blend two HEX colours

## Usage

``` r
blendColours(colour1, colour2, colour1Percentage = 0.5)
```

## Source

Code modified from <https://stackoverflow.com/questions/5560248>

## Arguments

- colour1:

  Character: HEX colour

- colour2:

  Character: HEX colour

- colour1Percentage:

  Character: percentage of colour 1 mixed in blended colour

## Value

Character representing an HEX colour

## Examples

``` r
psichomics:::blendColours("#3f83a3", "#f48000")
#> [1] "#9a8252"
```
