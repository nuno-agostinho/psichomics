# Parse an alternative splicing event from MISO

Parse an alternative splicing event from MISO

## Usage

``` r
parseMisoEvent(event)
```

## Arguments

- event:

  Data.frame containing only one event with at least 7 columns as
  retrieved from the alternative splicing annotation files from MISO
  (GFF3 files)

## Value

List with event attributes and junction positions for the exons (depends
on the events)

## Details

More information about MISO available at <http://miso.readthedocs.org>

## Examples

``` r
# example of alternative splicing event: skipped exon (SE)
event <- read.table(text = "
  chr1 SE gene 16854 18061 . - .
  chr1 SE mRNA 16854 18061 . - .
  chr1 SE exon 16854 17055 . - .
  chr1 SE exon 17233 17742 . - .
  chr1 SE exon 17915 18061 . - .
  chr1 SE mRNA 16854 18061 . - .
  chr1 SE exon 16854 17955 . - .
  chr1 SE exon 17915 18061 . - .")
psichomics:::parseMisoEvent(event)
#>   C1.start C1.end A1.start A1.end A2.start A2.end C2.start C2.end Program
#> 1    18061  17915    17742  17233       NA     NA    17055  16854    MISO
#>   Event.type Chromosome Strand
#> 1         SE       chr1      -
```
