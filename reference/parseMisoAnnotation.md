# Parse events from alternative splicing annotation

Parse events from alternative splicing annotation

## Usage

``` r
parseSuppaAnnotation(
  folder,
  types = c("SE", "AF", "AL", "MX", "A5", "A3", "RI"),
  genome = "hg19"
)

parseVastToolsAnnotation(
  folder,
  types = c("ALT3", "ALT5", "COMBI", "IR", "MERGE3m", "MIC", "EXSK", "MULTI"),
  genome = "Hsa",
  complexEvents = FALSE
)

parseMisoAnnotation(
  folder,
  types = c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI", "TandemUTR"),
  genome = "hg19"
)

parseMatsAnnotation(
  folder,
  types = c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI"),
  genome = "fromGTF",
  novelEvents = TRUE
)
```

## Arguments

- folder:

  Character: path to folder

- types:

  Character: type of events to retrieve (depends on the program of
  origin; see details)

- genome:

  Character: genome of interest (for instance, `hg19`; depends on the
  program of origin)

- complexEvents:

  Boolean: should complex events in A3SS and A5SS be parsed?

- novelEvents:

  Boolean: parse events detected due to novel splice sites

## Value

Retrieve data frame with events based on a given alternative splicing
annotation

## Details

Type of parsable events:

- Alternative 3' splice site

- Alternative 5' splice site

- Alternative first exon

- Alternative last exon

- Skipped exon (may include skipped micro-exons)

- Mutually exclusive exon

- Retained intron

- Tandem UTR

## See also

Other functions to prepare alternative splicing annotations:
[`prepareAnnotationFromEvents()`](https://nuno-agostinho.github.io/psichomics/reference/prepareAnnotationFromEvents.md)

## Examples

``` r
# Load sample files
folder <- "extdata/eventsAnnotSample/suppa_output/suppaEvents"
suppaOutput <- system.file(folder, package="psichomics")

suppa <- parseSuppaAnnotation(suppaOutput)
#> 
#> Retrieving SUPPA annotation...
#> 
#> Parsing SUPPA annotation...
# Load sample files
folder <- "extdata/eventsAnnotSample/VASTDB/Hsa/TEMPLATES"
vastToolsOutput <- system.file(folder, package="psichomics")

vast <- parseVastToolsAnnotation(vastToolsOutput)
#> 
#> Retrieving VAST-TOOLS annotation...
#> 
#> Parsing VAST-TOOLS annotation...
#> 
#> ALT3
#> 
#> ALT5
#> 
#> COMBI
#> 
#> EXSK
#> 
#> IR
#> 
#> MERGE3m
#> 
#> MIC
#> 
#> MULTI
# Load sample files
folder <- "extdata/eventsAnnotSample/miso_annotation"
misoOutput <- system.file(folder, package="psichomics")

miso <- parseMisoAnnotation(misoOutput)
#> 
#> Retrieving MISO annotation...
#> 
#> Parsing MISO annotation...
# Load sample files
folder <- "extdata/eventsAnnotSample/mats_output/ASEvents"
matsOutput <- system.file(folder, package="psichomics")

mats <- parseMatsAnnotation(matsOutput)
#> 
#> Retrieving rMATS annotation...
#> 
#> Parsing rMATS annotation...

# Do not parse novel events
mats <- parseMatsAnnotation(matsOutput, novelEvents=FALSE)
#> 
#> Retrieving rMATS annotation...
#> 
#> Parsing rMATS annotation...
```
