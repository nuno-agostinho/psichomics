# Render genetic information

Render genetic information

## Usage

``` r
renderGeneticInfo(
  output,
  info,
  species = NULL,
  assembly = NULL,
  grch37 = FALSE,
  eventDiagram = NULL,
  gene = NULL
)
```

## Arguments

- output:

  Shiny output

- info:

  Information as retrieved from Ensembl

- species:

  Character: species name

- assembly:

  Character: assembly version

- grch37:

  Boolean: use version GRCh37 of the genome?

- eventDiagram:

  Diagram of selected alternative splicing event

- ns:

  Namespace function

## Value

HTML elements to render gene, protein and transcript annotation
