# Parse junctions of an event from MISO according to event type

Parse junctions of an event from MISO according to event type

## Usage

``` r
parseMisoGeneric(event, validator, eventType, coord, plusIndex, minusIndex)

parseMisoSE(event)

parseMisoMXE(event)

parseMisoRI(event, strand)

parseMisoA5SS(event)

parseMisoA3SS(event, plusIndex, minusIndex)

parseMisoTandemUTR(event, minusIndex)

parseMisoAFE(event)

parseMisoALE(event)
```

## Arguments

- event:

  Data.frame containing only one event with at least 7 columns as
  retrieved from the alternative splicing annotation files from MISO
  (GFF3 files)

- validator:

  Character: valid elements for each event

- eventType:

  Character: event type (see details for available events)

- coord:

  Character: coordinate positions to fill

- plusIndex:

  Integer: index of the coordinates for a plus strand event

- minusIndex:

  Integer: index of the coordinates for a minus strand event

- strand:

  Character: positive-sense (`+`) or negative-sense `-` strand

## Value

List of parsed junctions

## Details

The following event types are available to be parsed:

- **SE** (exon skipping)

- **MXE** (mutually exclusive exon)

- **RI** (retained intron)

- **A5SS** (alternative 5' splice site)

- **A3SS** (alternative 3' splice site)

- **AFE** (alternative first exon)

- **ALE** (alternative last exon)

- **Tandem UTR**

## See also

[`parseMisoEvent()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoEvent.md)

## Examples

``` r
# skipped exon event (SE)
event <- read.table(text = "
  chr1 SE gene 16854 18061 . - .
  chr1 SE mRNA 16854 18061 . - .
  chr1 SE exon 16854 17055 . - .
  chr1 SE exon 17233 17742 . - .
  chr1 SE exon 17915 18061 . - .
  chr1 SE mRNA 16854 18061 . - .
  chr1 SE exon 16854 17955 . - .
  chr1 SE exon 17915 18061 . - .")
psichomics:::parseMisoSE(event)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end Program
#> 1    18061  17915    17742  17233       NA     NA    17055  16854    MISO
#>   Event.type Chromosome Strand
#> 1         SE       chr1      -

# mutually exclusive exon (MXE) event
event <- read.table(text = "
 chr1 MXE gene 764383 788090 . + .
 chr1 MXE mRNA 764383 788090 . + .
 chr1 MXE exon 764383 764484 . + .
 chr1 MXE exon 776580 776753 . + .
 chr1 MXE exon 787307 788090 . + .
 chr1 MXE mRNA 764383 788090 . + .
 chr1 MXE exon 764383 764484 . + .
 chr1 MXE exon 783034 783186 . + .
 chr1 MXE exon 787307 788090 . + .")
psichomics:::parseMisoMXE(event)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end Program
#> 1   764383 764484   776580 776753   783034 783186   787307 788090    MISO
#>   Event.type Chromosome Strand
#> 1        MXE       chr1      +

# retained intron (RI) event
event <- read.table(text = "
 chr1 RI gene 17233 17742 . - .
 chr1 RI mRNA 17233 17742 . - .
 chr1 RI exon 17233 17742 . - .
 chr1 RI mRNA 17233 17742 . - .
 chr1 RI exon 17233 17364 . - .
 chr1 RI exon 17601 17742 . - .")
psichomics:::parseMisoRI(event)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end Program
#> 1    17742  17601       NA     NA       NA     NA    17364  17233    MISO
#>   Event.type Chromosome Strand
#> 1         RI       chr1      -

# alternative 5' splice site (A5SS) event
event <- read.table(text = "
 chr1 A5SS gene 17233 17742 . - .
 chr1 A5SS mRNA 17233 17742 . - .
 chr1 A5SS exon 17233 17368 . - .
 chr1 A5SS exon 17526 17742 . - .
 chr1 A5SS mRNA 17233 17742 . - .
 chr1 A5SS exon 17233 17368 . - .
 chr1 A5SS exon 17606 17742 . - .")
psichomics:::parseMisoA5SS(event)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end Program
#> 1       NA     NA    17742  17526    17742  17606    17368  17233    MISO
#>   Event.type Chromosome Strand
#> 1       A5SS       chr1      -

# alternative 3' splice site (A3SS) event
event <- read.table(text = "
 chr1 A3SS gene 15796 16765 . - .
 chr1 A3SS mRNA 15796 16765 . - .
 chr1 A3SS exon 15796 15947 . - .
 chr1 A3SS exon 16607 16765 . - .
 chr1 A3SS mRNA 15796 16765 . - .
 chr1 A3SS exon 15796 15942 . - .
 chr1 A3SS exon 16607 16765 . - .")
psichomics:::parseMisoA3SS(event)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end Program
#> 1    16765  16607    15947  15796    15942  15796       NA     NA    MISO
#>   Event.type Chromosome Strand
#> 1       A3SS       chr1      -

# Tandem UTR event
event <- read.table(text = "
 chr19 TandemUTR gene  10663759  10664625  .  -  .
 chr19 TandemUTR mRNA  10663759  10664625  .  -  .
 chr19 TandemUTR exon  10663759  10664625  .  -  .
 chr19 TandemUTR mRNA  10664223  10664625  .  -  .
 chr19 TandemUTR exon  10664223  10664625  .  -  .")
psichomics:::parseMisoTandemUTR(event)
#>   C1.start C1.end A1.start   A1.end A2.start   A2.end C2.start C2.end Program
#> 1       NA     NA 10664625 10664223 10664625 10663759       NA     NA    MISO
#>   Event.type Chromosome Strand
#> 1  TandemUTR      chr19      -

# alternative first exon (AFE) event
event <- read.table(text = "
 chr12 AFE gene 57916659 57920171  .  +  .
 chr12 AFE mRNA 57919131 57920171  .  +  .
 chr12 AFE exon 57919131 57920171  .  +  .
 chr12 AFE mRNA 57916659 57918199  .  +  .
 chr12 AFE exon 57916659 57916794  .  +  .
 chr12 AFE exon 57917812 57917875  .  +  .
 chr12 AFE exon 57918063 57918199  .  +  .")
psichomics:::parseMisoAFE(event)
#>   C1.start C1.end A1.start   A1.end A2.start   A2.end C2.start C2.end Program
#> 1       NA     NA 57919131 57920171 57918063 57918199       NA     NA    MISO
#>   Event.type Chromosome Strand
#> 1        AFE      chr12      +

# alternative last exon (ALE) event
event <- read.table(text = "
 chr6 ALE gene 30620579 30822593  .  +  .
 chr6 ALE mRNA 30822190 30822593  .  +  .
 chr6 ALE exon 30822190 30822593  .  +  .
 chr6 ALE mRNA 30620579 30620982  .  +  .
 chr6 ALE exon 30620579 30620982  .  +  .")
psichomics:::parseMisoALE(event)
#>   C1.start C1.end A1.start   A1.end A2.start   A2.end C2.start C2.end Program
#> 1       NA     NA 30822190 30822593 30620579 30620982       NA     NA    MISO
#>   Event.type Chromosome Strand
#> 1        ALE       chr6      +
```
