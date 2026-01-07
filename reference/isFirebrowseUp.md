# Check if [FireBrowse API](http://firebrowse.org/api-docs/) is running

Check if [FireBrowse API](http://firebrowse.org/api-docs/) is running

## Usage

``` r
isFirebrowseUp()
```

## Value

Invisible `TRUE` if the [FireBrowse
API](http://firebrowse.org/api-docs/) is working; otherwise, raises a
warning with the status code and a brief explanation.

## See also

Other functions associated with TCGA data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`getTCGAdataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md),
[`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md),
[`parseTCGAsampleTypes()`](https://nuno-agostinho.github.io/psichomics/reference/parseTCGAsampleInfo.md)

## Examples

``` r
isFirebrowseUp()
#> [1] TRUE
```
