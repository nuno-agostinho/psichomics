# Load alternative splicing annotation from `AnnotationHub`

Load alternative splicing annotation from `AnnotationHub`

## Usage

``` r
loadAnnotation(annotation, cache = getAnnotationHubOption("CACHE"))
```

## Arguments

- annotation:

  Character: annotation to load

- cache:

  Character: path to `AnnotationHub` cache (used to load alternative
  splicing event annotation)

## Value

List of data frames containing the alternative splicing annotation per
event type

## See also

Other functions for PSI quantification:
[`filterPSI()`](https://nuno-agostinho.github.io/psichomics/reference/filterPSI.md),
[`getSplicingEventTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventTypes.md),
[`listSplicingAnnotations()`](https://nuno-agostinho.github.io/psichomics/reference/listSplicingAnnotations.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md),
[`quantifySplicing()`](https://nuno-agostinho.github.io/psichomics/reference/quantifySplicing.md)

## Examples

``` r
human <- listSplicingAnnotations(species="Homo sapiens")[[1]]
if (FALSE) { # \dontrun{
annot <- loadAnnotation(human)
} # }
```
