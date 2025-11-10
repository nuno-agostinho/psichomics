# Calculate inclusion levels using alternative splicing event annotation and junction quantification for many samples

Calculate inclusion levels using alternative splicing event annotation
and junction quantification for many samples

## Usage

``` r
calculateInclusionLevels(
  eventType,
  junctionQuant,
  annotation,
  minReads = 10,
  onlyReturnASeventNames = FALSE
)
```

## Arguments

- eventType:

  Character: type of the alternative event to calculate

- junctionQuant:

  Matrix: junction quantification with samples as columns and junctions
  as rows

- annotation:

  Data.frame: alternative splicing annotation related to event type

- minReads:

  Integer: minimum of total reads required to consider the
  quantification as valid

## Value

Matrix with inclusion levels
