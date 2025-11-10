# Parse junctions of an event from SUPPA

Parse junctions of an event from SUPPA

## Usage

``` r
parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)

parseSuppaSE(junctions, strand)

parseSuppaRI(junctions, strand)

parseSuppaALE(junctions, strand)

parseSuppaAFE(junctions, strand)

parseSuppaMXE(junctions, strand)

parseSuppaA3SS(junctions, strand)

parseSuppaA5SS(junctions, strand)
```

## Arguments

- junctions:

  List of integers: exon-exon junctions of an event

- strand:

  Character: positive-sense (`+`) or negative-sense (`-`) strand

- coords:

  Character: coordinate positions to fill

- plus_pos:

  Integer: index of the coordinates for a plus strand event

- minus_pos:

  Integer: index of the coordinates for a minus strand event

## Value

Data frame of parsed junctions

## Details

The following event types are available to be parsed:

- **SE** (exon skipping)

- **RI** (retained intron)

- **MXE** (mutually exclusive exons)

- **A5SS** (alternative 5' splice site)

- **A3SS** (alternative 3' splice site)

- **ALE** (alternative last exon)

- **AFE** (alternative first exon)

## See also

[`parseSuppaEvent()`](https://nuno-agostinho.github.io/psichomics/reference/parseSuppaEvent.md)

## Examples

``` r
# Parse generic event (in this case, an exon skipping event)
junctions <- read.table(text = "169768099 169770024 169770112 169771762")
coords <- c("C1.end", "A1.start", "A1.end", "C2.start")
plus  <- 1:4
minus <- 1:4
psichomics:::parseSuppaGeneric(junctions, strand = "+", coords, plus, minus)
#>   C1.start    C1.end  A1.start    A1.end A2.start A2.end  C2.start C2.end
#> 1       NA 169768099 169770024 169770112       NA     NA 169771762     NA

junctions <- read.table(text = "169768099 169770024 169770112 169771762")
psichomics:::parseSuppaSE(junctions, "+")
#>   C1.start    C1.end  A1.start    A1.end A2.start A2.end  C2.start C2.end
#> 1       NA 169768099 169770024 169770112       NA     NA 169771762     NA

junctions <- read.table(text = "196709749 196709922 196711005 196711181")
psichomics:::parseSuppaRI(junctions, "+")
#>    C1.start    C1.end A1.start A1.end A2.start A2.end  C2.start    C2.end
#> 1 196709749 196709922       NA     NA       NA     NA 196711005 196711181

junctions <- read.table(
    text = "24790610 24792494 24792800 24790610 24795476 24795797")
psichomics:::parseSuppaALE(junctions, "+")
#>   C1.start   C1.end A1.start   A1.end A2.start   A2.end C2.start C2.end
#> 1       NA 24790610 24792494 24792800 24795476 24795797       NA     NA

junctions <- read.table(
    text = "169763871 169764046 169767998 169764550 169765124 169767998")
psichomics:::parseSuppaAFE(junctions, "+")
#>   C1.start C1.end  A1.start    A1.end  A2.start    A2.end  C2.start C2.end
#> 1       NA     NA 169763871 169764046 169764550 169765124 169767998     NA

junctions <- read.table(
    text = "202060671 202068453 202068489 202073793 202060671 202072798 202072906 202073793")
psichomics:::parseSuppaMXE(junctions, "+")
#>   C1.start    C1.end  A1.start    A1.end  A2.start    A2.end  C2.start C2.end
#> 1       NA 202060671 202068453 202068489 202072798 202072906 202073793     NA

junctions <- read.table(text = "169772450 169773216 169772450 169773253")
psichomics:::parseSuppaA3SS(junctions, "+")
#>   C1.start    C1.end  A1.start A1.end  A2.start A2.end C2.start C2.end
#> 1       NA 169772450 169773216     NA 169773253     NA       NA     NA

junctions <- read.table(text = "50193276 50197008 50192997 50197008")
psichomics:::parseSuppaA5SS(junctions, "+")
#>   C1.start C1.end A1.start   A1.end A2.start   A2.end C2.start C2.end
#> 1       NA     NA       NA 50193276       NA 50192997 50197008     NA
```
