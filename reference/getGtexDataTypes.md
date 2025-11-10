# Get GTEx data information

Get GTEx data information

## Usage

``` r
getGtexDataTypes()

getGtexReleases()
```

## Value

GTEx data information

## See also

Other functions associated with GTEx data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`getGtexTissues()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexTissues.md),
[`loadGtexData()`](https://nuno-agostinho.github.io/psichomics/reference/loadGtexData.md)

## Examples

``` r
getGtexDataTypes()
#>       Sample attributes      Subject phenotypes         Gene expression 
#>            "sampleInfo"           "subjectInfo"              "geneExpr" 
#> Junction quantification 
#>         "junctionQuant" 
getGtexReleases()
#> [1] 8 7 6 4
```
