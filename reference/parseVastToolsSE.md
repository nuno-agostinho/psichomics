# Parse junctions of an event from VAST-TOOLS according to event type

Parse junctions of an event from VAST-TOOLS according to event type

## Usage

``` r
parseVastToolsSE(junctions)

parseVastToolsRI(junctions, strand)

parseVastToolsA3SS(junctions)

parseVastToolsA5SS(junctions)
```

## Arguments

- junctions:

  Data.frame or matrix: exon-exon junctions of alternative splicing
  events (it must have 4 columns)

- strand:

  Character: positive (+) or negative (-) strand

## Value

List of parsed junctions

## Details

The following event types are available to be parsed:

- **SE** (skipped exon)

- **RI** (retained intron)

- **A5SS** (alternative 5' splice site)

- **A3SS** (alternative 3' splice site)

## See also

[`parseVastToolsEvent()`](https://nuno-agostinho.github.io/psichomics/reference/parseVastToolsEvent.md)

## Examples

``` r
junctions <- read.table(text = "41040823 41046768 41046903 41051785")
psichomics:::parseVastToolsSE(junctions)
#>   C1.start   C1.end A1.start   A1.end A2.start A2.end C2.start C2.end Strand
#> 1       NA 41040823 41046768 41046903       NA     NA 41051785     NA      +

# these functions are vectorised!
junctions <- read.table(text = "41040823 41046768 41046903 41051785
                                58864658 58864693 58864294 58864563")
psichomics:::parseVastToolsSE(junctions)
#>   C1.start   C1.end A1.start   A1.end A2.start A2.end C2.start C2.end Strand
#> 1       NA 41040823 41046768 41046903       NA     NA 41051785     NA      +
#> 2       NA 58864658 58864294 58864693       NA     NA 58864563     NA      -

junctions <- read.table(text = "58864658 58864693 58864294 58864563")
psichomics:::parseVastToolsRI(junctions, strand = "+")
#>   C1.start   C1.end A1.start A1.end A2.start A2.end C2.start   C2.end Strand
#> 1 58864658 58864693       NA     NA       NA     NA 58864294 58864563      +

junctions <- rbind(
    c(36276385, list(c(36277798, 36277315)), 36277974),
    c(7133604, 7133377, list(c(7133474, 7133456)))
)
psichomics:::parseVastToolsA3SS(junctions)
#>   C1.start   C1.end A1.start A1.end A2.start   A2.end C2.start C2.end Strand
#> 1       NA 36276385 36277798     NA 36277315 36277974       NA     NA      +
#> 2       NA  7133604  7133474     NA  7133456  7133377       NA     NA      -

junctions <- rbind(
    c(74650610, list(c(74650654, 74650658)), 74650982),
    c(list(c(49557666, 49557642), 49557746, 49557470))
)
psichomics:::parseVastToolsA5SS(junctions)
#>   C1.start C1.end A1.start   A1.end A2.start   A2.end C2.start C2.end Strand
#> 1       NA     NA       NA 74650654 74650610 74650658 74650982     NA      +
#> 2       NA     NA       NA 49557666 49557746 49557642 49557470     NA      -
```
