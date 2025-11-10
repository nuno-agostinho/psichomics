# Parses an alternative splicing event from VAST-TOOLS

Parses an alternative splicing event from VAST-TOOLS

## Usage

``` r
parseVastToolsEvent(event)
```

## Arguments

- event:

  Data.frame: VAST-TOOLS event containing gene symbol, event ID, length,
  junctions coordinates, event type and inclusion levels for both
  samples

## Value

List with the event attributes (chromosome, strand, event type and the
position of the exon boundaries)

## Details

Junctions are parsed from

## Note

Only supports to parse one event at a time.

## Examples

``` r
event <- read.table(text =
"NFYA HsaEX0042823 chr6:41046768-41046903 136 chr6:41040823,41046768-41046903,41051785 C2 0 N 0 N"
)
psichomics:::parseVastToolsEvent(event)
#>      Program Gene.symbol     Event.ID Event.type Chromosome Inclusion.level.A
#> 1 VAST-TOOLS        NFYA HsaEX0042823         SE       chr6                 0
#>   Inclusion.level.B C1.start   C1.end A1.start   A1.end A2.start A2.end
#> 1                 0       NA 41040823 41046768 41046903       NA     NA
#>   C2.start C2.end Strand
#> 1 41051785     NA      +
```
