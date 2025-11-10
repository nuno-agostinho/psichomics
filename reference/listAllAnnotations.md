# List alternative splicing annotation files available, as well as custom annotation

List alternative splicing annotation files available, as well as custom
annotation

## Usage

``` r
listAllAnnotations(...)
```

## Arguments

- ...:

  Custom annotation loaded

## Value

Named character vector with splicing annotation files available

## Examples

``` r
psichomics:::listAllAnnotations()
#> $`MISO, VAST-TOOLS, UCSC`
#> Human hg19 (2016-10-11) Human hg19 (2017-10-20) Human hg38 (2018-04-30) 
#>               "AH51461"               "AH60272"               "AH63657" 
#> 
#> $`VAST-TOOLS`
#>                            Human hg19 (2021-06-15) 
#>                                          "AH95569" 
#>                            Human hg38 (2021-06-15) 
#>                                          "AH95570" 
#>                      Mus musculus mm9 (2021-06-15) 
#>                                          "AH95571" 
#>                     Mus musculus mm10 (2021-06-15) 
#>                                          "AH95572" 
#>                    Bos taurus bosTau6 (2021-06-15) 
#>                                          "AH95573" 
#>                 Gallus gallus galGal3 (2021-06-15) 
#>                                          "AH95574" 
#>                 Gallus gallus galGal4 (2021-06-15) 
#>                                          "AH95575" 
#>            Xenopus tropicalis xenTro3 (2021-06-15) 
#>                                          "AH95576" 
#>                  Danio rerio danRer10 (2021-06-15) 
#>                                          "AH95577" 
#>     Branchiostoma lanceolatum braLan2 (2021-06-15) 
#>                                          "AH95578" 
#> Strongylocentrotus purpuratus strPur4 (2021-06-15) 
#>                                          "AH95579" 
#>           Drosophila melanogaster dm6 (2021-06-15) 
#>                                          "AH95580" 
#>            Strigamia maritima strMar1 (2021-06-15) 
#>                                          "AH95581" 
#>           Caenorhabditis elegans ce11 (2021-06-15) 
#>                                          "AH95582" 
#>       Schmidtea mediterranea schMed31 (2021-06-15) 
#>                                          "AH95583" 
#>        Nematostella vectensis nemVec1 (2021-06-15) 
#>                                          "AH95584" 
#>         Arabidopsis thaliana araTha10 (2021-06-15) 
#>                                          "AH95585" 
#> 
#> $`Custom annotation`
#> Load annotation from file... 
#>             "loadAnnotation" 
#> 
```
