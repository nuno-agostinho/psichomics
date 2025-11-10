# Returns the date format used by the FireBrowse API

Returns the date format used by the FireBrowse API

## Usage

``` r
getFirebrowseDateFormat()
```

## Value

Named list with date formats from FireBrowse API

## Examples

``` r
format <- psichomics:::getFirebrowseDateFormat()

# date format to use in a query to FireBrowse API
format$query
#> [1] "%Y_%m_%d"

# date format to parse a date in a response from FireBrowse API
format$response
#> [1] "%d %b %Y"
```
