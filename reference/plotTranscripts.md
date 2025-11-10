# Plot transcripts

Plot transcripts

## Usage

``` r
plotTranscripts(
  info,
  eventPosition = NULL,
  event = NULL,
  eventData = NULL,
  shiny = FALSE
)
```

## Arguments

- info:

  Information retrieved from Ensembl

- eventPosition:

  Numeric: coordinates of the alternative splicing event (ignored if
  `event` is set)

- event:

  Character: identifier of the alternative splicing event to plot

- eventData:

  Object containing event information to be parsed

- shiny:

  Boolean: is the function running in a Shiny session?

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## See also

Other functions to retrieve external information:
[`ensemblToUniprot()`](https://nuno-agostinho.github.io/psichomics/reference/ensemblToUniprot.md),
[`plotProtein()`](https://nuno-agostinho.github.io/psichomics/reference/plotProtein.md),
[`queryEnsemblByGene()`](https://nuno-agostinho.github.io/psichomics/reference/queryEnsemblByGene.md)

## Examples

``` r
event <- "SE_12_-_7985318_7984360_7984200_7982602_SLC2A14"
info  <- queryEnsemblByEvent(event, species="human", assembly="hg19")
if (FALSE) { # \dontrun{
plotTranscripts(info, event=event)
} # }
```
