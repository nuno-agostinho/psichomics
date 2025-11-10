# Parse file according to its format

Parse file according to its format

## Usage

``` r
parseFile(format, file, ..., verbose = FALSE)
```

## Arguments

- format:

  Environment: format of the file

- file:

  Character: file to load

- ...:

  Extra parameters passed to
  [fread](https://rdatatable.gitlab.io/data.table/reference/fread.html)

- verbose:

  Boolean: detail step while parsing?

## Value

Data frame with the loaded file

## Details

The resulting data frame includes the attribute `tablename` with the
name of the data frame
