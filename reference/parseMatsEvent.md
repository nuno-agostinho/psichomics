# Parse alternative splicing events from MATS

Parse alternative splicing events from MATS

## Usage

``` r
parseMatsEvent(event, event_type)
```

## Arguments

- event:

  Data frame row: MATS splicing event

- event_type:

  Character: Type of event to parse (see details)

## Value

List containing the event attributes and junctions

## Details

The following event types can be parsed:

- **SE**: Skipped exon

- **MXE**: Mutually exclusive exons

- **RI**:Retained intron

- **A3SS**: Alternative 3' splice site

- **A5SS**: Alternative 5' splice site

## Examples

``` r
# MATS event (alternative 3' splice site)
event <- read.table(text = "
     2 ENSG00000166012 TAF1D chr11 - 93466515 93466671 93466515 93466563 93467790 93467826
     5 ENSG00000166012 TAF1D chr11 - 93466515 93466671 93466515 93466585 93467790 93467826
     6 ENSG00000166012 TAF1D chr11 - 93466515 93466585 93466515 93466563 93467790 93467826
")
psichomics:::parseMatsEvent(event, "A3SS")
#>   Program  Gene Chromosome Strand Event.type      Event.ID C1.start   C1.end
#> 1    MATS TAF1D      chr11      -       A3SS unassigned_id 93467826 93467790
#> 2    MATS TAF1D      chr11      -       A3SS unassigned_id 93467826 93467790
#> 3    MATS TAF1D      chr11      -       A3SS unassigned_id 93467826 93467790
#>   A1.start A1.end A2.start   A2.end C2.start C2.end
#> 1 93466671     NA 93466563 93466515       NA     NA
#> 2 93466671     NA 93466585 93466515       NA     NA
#> 3 93466585     NA 93466563 93466515       NA     NA
```
