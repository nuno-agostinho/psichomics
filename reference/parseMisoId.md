# Parse MISO's alternative splicing event identifier

Parse MISO's alternative splicing event identifier

## Usage

``` r
parseMisoId(id)
```

## Arguments

- id:

  Character: MISO alternative splicing event identifier

## Value

Character with the parsed ID

## Examples

``` r
id <- paste0(
    "ID=ENSMUSG00000026150.chr1:82723803:82723911:+@chr1:82724642:82724813:",
    "+@chr1:82725791:82726011:+.B;Parent=ENSMUSG00000026150.chr1:82723803:",
    "82723911:+@chr1:82724642:82724813:+@chr1:82725791:82726011:+")
psichomics:::parseMisoId(id)
#> [1] "ENSMUSG00000026150.chr1:82723803:82723911:+@chr1:82724642:82724813:+@chr1:82725791:82726011:+.B"
```
