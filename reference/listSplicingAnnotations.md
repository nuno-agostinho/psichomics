# List alternative splicing annotations

List alternative splicing annotations

## Usage

``` r
listSplicingAnnotations(
  species = NULL,
  assembly = NULL,
  date = NULL,
  cache = getAnnotationHubOption("CACHE"),
  group = FALSE
)
```

## Arguments

- species:

  Character: filter results by species (regular expression)

- assembly:

  Character: filter results by assembly (regular expression)

- date:

  Character: filter results by date (regular expression)

- cache:

  Character: path to `AnnotationHub` cache (used to load alternative
  splicing event annotation)

- group:

  Boolean: group values based on data provider?

## Value

Named character vector with splicing annotation names

## See also

Other functions for PSI quantification:
[`filterPSI()`](https://nuno-agostinho.github.io/psichomics/reference/filterPSI.md),
[`getSplicingEventTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventTypes.md),
[`loadAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/loadAnnotation.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md),
[`quantifySplicing()`](https://nuno-agostinho.github.io/psichomics/reference/quantifySplicing.md)

## Examples

``` r
listSplicingAnnotations() # Return all alternative splicing annotations
#>                                            Human hg19 (2016-10-11) 
#>                                                          "AH51461" 
#>                                            Human hg19 (2017-10-20) 
#>                                                          "AH60272" 
#>                                            Human hg38 (2018-04-30) 
#>                                                          "AH63657" 
#>                            Human hg19 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95569" 
#>                            Human hg38 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95570" 
#>                      Mus musculus mm9 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95571" 
#>                     Mus musculus mm10 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95572" 
#>                    Bos taurus bosTau6 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95573" 
#>                 Gallus gallus galGal3 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95574" 
#>                 Gallus gallus galGal4 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95575" 
#>            Xenopus tropicalis xenTro3 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95576" 
#>                  Danio rerio danRer10 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95577" 
#>     Branchiostoma lanceolatum braLan2 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95578" 
#> Strongylocentrotus purpuratus strPur4 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95579" 
#>           Drosophila melanogaster dm6 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95580" 
#>            Strigamia maritima strMar1 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95581" 
#>           Caenorhabditis elegans ce11 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95582" 
#>       Schmidtea mediterranea schMed31 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95583" 
#>        Nematostella vectensis nemVec1 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95584" 
#>         Arabidopsis thaliana araTha10 from VAST-TOOLS (2021-06-15) 
#>                                                          "AH95585" 
listSplicingAnnotations(assembly="hg19") # Search for hg19 annotation
#>                 Human hg19 (2016-10-11)                 Human hg19 (2017-10-20) 
#>                               "AH51461"                               "AH60272" 
#> Human hg19 from VAST-TOOLS (2021-06-15) 
#>                               "AH95569" 
listSplicingAnnotations(assembly="hg38") # Search for hg38 annotation
#>                 Human hg38 (2018-04-30) Human hg38 from VAST-TOOLS (2021-06-15) 
#>                               "AH63657"                               "AH95570" 
listSplicingAnnotations(date="201(7|8)") # Search for 2017 or 2018 annotation
#> Human hg19 (2017-10-20) Human hg38 (2018-04-30) 
#>               "AH60272"               "AH63657" 
```
