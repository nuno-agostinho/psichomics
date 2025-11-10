# Parse junctions of an alternative splicing event from MATS according to event type

Parse junctions of an alternative splicing event from MATS according to
event type

## Usage

``` r
parseMatsGeneric(junctions, strand, coords, plus_pos, minus_pos)

parseMatsSE(junctions, strand)

parseMatsMXE(junctions, strand)

parseMatsRI(junctions, strand)

parseMatsA3SS(junctions, strand)

parseMatsA5SS(junctions, strand)

parseMatsAFE(junctions, strand)

parseMatsALE(junctions, strand)
```

## Arguments

- junctions:

  Integer: event's junctions

- strand:

  Character: strand of the event

- coords:

  Character: names of the alternative splicing coordinates

- plus_pos:

  Integer: match of each junction in the respective coordinate for the
  plus strand

- minus_pos:

  Integer: match of each junction in the respective coordinate for the
  minus strand

## Value

Data frame with parsed junctions

## Details

The following event types are ready to be parsed:

- **SE** (skipped exon)

- **MXE** (mutually exclusive exon)

- **RI** (retained intron)

- **A5SS** (alternative 5' splice site)

- **A3SS** (alternative 3' splice site)

- **AFE** (alternative first exon)

- **ALE** (alternative last exon)

You can use `parseMatsGeneric` to parse other event types.

## See also

[`parseMatsEvent()`](https://nuno-agostinho.github.io/psichomics/reference/parseMatsEvent.md)

## Examples

``` r
# Parse generic event (in this case, an exon skipping event)
junctions <- read.table(text=
    "79685787 79685910 79685796 79685910 79679566 79679751")
coords <- c("A1.start", "A1.end",
            "C1.start", "C1.end",
            "C2.start", "C2.end")
plus  <- c(1:6)
minus <- c(2:1, 6:3)
psichomics:::parseMatsGeneric(junctions, strand = "+", coords, plus, minus)
#>   C1.start   C1.end A1.start   A1.end A2.start A2.end C2.start   C2.end
#> 1 79685796 79685910 79685787 79685910       NA     NA 79679566 79679751

# Parse exon skipping event
junctions <- read.table(text=
    "79685787 79685910 79685796 79685910 79679566 79679751")
psichomics:::parseMatsSE(junctions, strand = "+")
#>   C1.start   C1.end A1.start   A1.end A2.start A2.end C2.start   C2.end
#> 1 79685796 79685910 79685787 79685910       NA     NA 79679566 79679751

# Parse mutually exclusive exon event
junctions <- read.table(text=
"158282161 158282276 158282689 158282804 158281047 158281295 158283950 158284199")
psichomics:::parseMatsMXE(junctions, strand = "+")
#>    C1.start    C1.end  A1.start    A1.end  A2.start    A2.end  C2.start
#> 1 158281047 158281295 158282161 158282276 158282689 158282804 158283950
#>      C2.end
#> 1 158284199

# Parse retained intron event
junctions <- read.table(text=
    "15929853 15932100 15929853 15930016 15930687 15932100")
psichomics:::parseMatsRI(junctions, strand = "+")
#>   C1.start   C1.end A1.start A1.end A2.start A2.end C2.start   C2.end
#> 1 15929853 15930016       NA     NA       NA     NA 15930687 15932100

# Parse alternative 3' splicing site event
junctions <- read.table(text=
    "79685787 79685910 79685796 79685910 79679566 79679751")
psichomics:::parseMatsA3SS(junctions, strand = "+")
#>   C1.start   C1.end A1.start A1.end A2.start   A2.end C2.start C2.end
#> 1 79679566 79679751 79685787     NA 79685796 79685910       NA     NA

# Parse alternative 5' splicing site event
junctions <- read.table(text=
    "102884421 102884501 102884421 102884489 102884812 102885881")
psichomics:::parseMatsA5SS(junctions, strand = "+")
#>   C1.start C1.end A1.start    A1.end  A2.start    A2.end  C2.start    C2.end
#> 1       NA     NA       NA 102884501 102884421 102884489 102884812 102885881

# Parse alternative first exon event
junctions <- read.table(text=
    "16308723 16308879 16308967 16309119 16314269 16314426")
psichomics:::parseMatsAFE(junctions, strand = "+")
#>   C1.start C1.end A1.start   A1.end A2.start   A2.end C2.start   C2.end
#> 1       NA     NA 16308723 16308879 16308967 16309119 16314269 16314426

# Parse alternative last exon event
junctions <- read.table(text=
    "111858645 111858828 111851063 111851921 111850441 111850543")
psichomics:::parseMatsAFE(junctions, strand = "+")
#>   C1.start C1.end  A1.start    A1.end  A2.start    A2.end  C2.start    C2.end
#> 1       NA     NA 111858645 111858828 111851063 111851921 111850441 111850543
```
