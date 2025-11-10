# Parses splicing events of a specific event type from SUPPA

Parses splicing events of a specific event type from SUPPA

## Usage

``` r
parseSuppaEvent(event)
```

## Arguments

- event:

  Character vector: Splicing event attributes and junction positions

## Value

List with the event attributes (chromosome, strand, event type and the
position of the exon boundaries)

## Details

More information about SUPPA available at
<https://bitbucket.org/regulatorygenomicsupf/suppa>

The following event types are available to be parsed:

- **SE** (skipped exon)

- **RI** (retained intron)

- **MX** (mutually exclusive exons)

- **A5** (alternative 5' splice site)

- **A3** (alternative 3' splice site)

- **AL** (alternative last exon)

- **AF** (alternative first exon)

## Note

It only allows to parse one event type at once.

## Examples

``` r
event <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
psichomics:::parseSuppaEvent(event)
#>   Program            Gene
#> 1   SUPPA ENSG00000000419
#>                                                      Event.ID Chromosome
#> 1 ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-         20
#>   Event.type Strand C1.start   C1.end A1.start A1.end A2.start A2.end C2.start
#> 1       A3SS      -       NA 49557642 49557470     NA 49557492     NA       NA
#>   C2.end
#> 1     NA
```
