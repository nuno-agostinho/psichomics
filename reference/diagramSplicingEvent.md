# Prepare SVG diagram of alternative splicing events

Prepare SVG diagram of alternative splicing events

## Usage

``` r
diagramSplicingEvent(
  parsed,
  type,
  class = "pull-right",
  style = NULL,
  showText = TRUE,
  showPath = TRUE,
  showAlternative1 = TRUE,
  showAlternative2 = TRUE,
  constitutiveWidth = NULL,
  alternativeWidth = NULL,
  intronWidth = NULL,
  constitutiveFill = "lightgray",
  constitutiveStroke = "darkgray",
  alternative1Fill = "#ffb153",
  alternative1Stroke = "#faa000",
  alternative2Fill = "#caa06c",
  alternative2Stroke = "#9d7039"
)
```

## Arguments

- parsed:

  Alternative splicing event

- type:

  Character: alternative splicing event type

- class:

  Character: class of SVG parent tag

- style:

  Character: style of SVG parent tag

- showText:

  Boolean: display coordinates and length (if available)

- showPath:

  Boolean: display alternative splicing junctions

- showAlternative1:

  Boolean: show alternative exon 1 and respective splicing junctions and
  text?

- showAlternative2:

  Boolean: show alternative exon 2 and respective splicing junctions and
  text? (only related with mutually exclusive exons)

- constitutiveWidth:

  Numeric: width of constitutive exon(s)

- alternativeWidth:

  Numeric: width of alternative exon(s)

- intronWidth:

  Numeric: width of intron's representation

- constitutiveFill:

  Character: fill colour of constitutive exons

- constitutiveStroke:

  Character: stroke colour of constitutive exons

- alternative1Fill:

  Character: fill colour of alternative exon 1

- alternative1Stroke:

  Character: stroke colour of alternative exon 1

- alternative2Fill:

  Character: fill colour of alternative exon 2

- alternative2Stroke:

  Character: stroke colour of alternative exon 2

## Value

Diagrams per alternative splicing event in SVG
