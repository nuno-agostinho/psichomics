# Get GTEx tissues from given GTEx sample attributes

Get GTEx tissues from given GTEx sample attributes

## Usage

``` r
getGtexTissues(folder = getDownloadsFolder(), release = getGtexReleases()[[1]])
```

## Arguments

- folder:

  Character: folder containing data

- release:

  Numeric: GTEx data release to load

## Value

Character: available tissues

## See also

Other functions associated with GTEx data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`getGtexDataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexDataTypes.md),
[`loadGtexData()`](https://nuno-agostinho.github.io/psichomics/reference/loadGtexData.md)

## Examples

``` r
if (FALSE) { # \dontrun{
getGtexTissues()
} # }
```
