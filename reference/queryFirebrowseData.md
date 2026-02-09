# Query the FireBrowse API for TCGA data

Query the FireBrowse API for TCGA data

## Usage

``` r
queryFirebrowseData(
  format = "json",
  date = NULL,
  cohort = NULL,
  data_type = NULL,
  tool = NULL,
  platform = NULL,
  center = NULL,
  level = NULL,
  protocol = NULL,
  page = NULL,
  page_size = NULL,
  sort_by = NULL
)
```

## Arguments

- format:

  Character: response format as `JSON`, `CSV` or `TSV`

- date:

  Character: dates of the data retrieval by FireBrowse (by default, it
  uses the most recent data available)

- cohort:

  Character: abbreviation of the cohorts (by default, returns data for
  all cohorts)

- data_type:

  Character: data types (optional)

- tool:

  Character: data produced by the selected FireBrowse tools (optional)

- platform:

  Character: data generation platforms (optional)

- center:

  Character: data generation centres (optional)

- level:

  Integer: data levels (optional)

- protocol:

  Character: sample characterization protocols (optional)

- page:

  Integer: page of the results to return (optional)

- page_size:

  Integer: number of records per page of results (optional)

- sort_by:

  String: column used to sort the data (by default, sort by cohort)

## Value

Response from the FireBrowse API (it needs to be parsed)

## Examples

``` r
cohort <- getTCGAcohorts()[1]
psichomics:::queryFirebrowseData(cohort = names(cohort),
                                 data_type = "mRNASeq")
#> Response [http://firebrowse.org/api/v1/Archives/StandardData?format=json&date=2016_01_28&cohort=ACC&data_type=mRNASeq]
#>   Date: 2026-02-09 08:03
#>   Status: 200
#>   Content-Type: application/json
#>   Size: 12.1 kB
#> {
#>   "StandardData": [
#>     {
#>       "active": true,
#>       "center": "unc.edu",
#>       "cohort": "ACC",
#>       "data_type": "mRNASeq",
#>       "date": "Thu, 28 Jan 2016 00:00:00 GMT",
#>       "filesize": 16560553,
#>       "level": 3,
#> ...

# Querying for data from a specific date
dates <- getTCGAdates()
dates <- format(dates, psichomics:::getFirebrowseDateFormat()$query)

psichomics:::queryFirebrowseData(date = dates[2], cohort = names(cohort))
#> Response [http://firebrowse.org/api/v1/Archives/StandardData?format=json&date=2015_11_01&cohort=ACC]
#>   Date: 2026-02-09 08:03
#>   Status: 200
#>   Content-Type: application/json
#>   Size: 37.2 kB
#> {
#>   "StandardData": [
#>     {
#>       "active": true,
#>       "center": "unc.edu",
#>       "cohort": "ACC",
#>       "data_type": "mRNASeq",
#>       "date": "Sun, 01 Nov 2015 00:00:00 GMT",
#>       "level": 3,
#>       "platform": "illuminahiseq_rnaseqv2",
#> ...
```
