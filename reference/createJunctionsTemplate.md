# Creates a template of alternative splicing junctions

Creates a template of alternative splicing junctions

## Usage

``` r
createJunctionsTemplate(
  nrow,
  program = character(0),
  event.type = character(0),
  chromosome = character(0),
  strand = character(0),
  id = character(0)
)
```

## Arguments

- nrow:

  Integer: row number

- program:

  Character: program used to get the junctions

- event.type:

  Character: event type

- chromosome:

  Character: chromosome

- strand:

  Character: positive-sense (`+`) or negative-sense (`-`) strand

- id:

  Character: event identifiers

## Value

A data frame with the junctions coordinate names pre-filled with `NA`

## Examples

``` r
psichomics:::createJunctionsTemplate(nrow = 8)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end
#> 1       NA     NA       NA     NA       NA     NA       NA     NA
#> 2       NA     NA       NA     NA       NA     NA       NA     NA
#> 3       NA     NA       NA     NA       NA     NA       NA     NA
#> 4       NA     NA       NA     NA       NA     NA       NA     NA
#> 5       NA     NA       NA     NA       NA     NA       NA     NA
#> 6       NA     NA       NA     NA       NA     NA       NA     NA
#> 7       NA     NA       NA     NA       NA     NA       NA     NA
#> 8       NA     NA       NA     NA       NA     NA       NA     NA
```
