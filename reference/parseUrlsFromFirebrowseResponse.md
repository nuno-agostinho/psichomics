# Retrieve URLs from a response to a FireBrowse data query

Retrieve URLs from a response to a FireBrowse data query

## Usage

``` r
parseUrlsFromFirebrowseResponse(res)
```

## Arguments

- res:

  Response from [`httr::GET`](https://httr.r-lib.org/reference/GET.html)
  to a FireBrowse data query

## Value

Named character with URLs

## Examples

``` r
res <- psichomics:::queryFirebrowseData(cohort = "ACC")
url <- psichomics:::parseUrlsFromFirebrowseResponse(res)
```
