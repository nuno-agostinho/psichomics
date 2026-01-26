# Load file based on its format

Tries to recognise the file format and parses the content of the given
file accordingly.

## Usage

``` r
loadFile(
  file,
  formats = loadFileFormats(),
  ...,
  verbose = FALSE,
  multiple = FALSE
)
```

## Arguments

- file:

  Character: file to parse

- formats:

  List of file formats to check

- ...:

  Extra parameters passed to
  [fread](https://rdrr.io/pkg/data.table/man/fread.html)

- verbose:

  Boolean: detail steps while parsing

- multiple:

  Boolean: expect more than one file?

## Value

Data frame with the contents of the given file if the file format is
recognised; otherwise, returns `NULL`

## Details

The resulting data frame includes the attribute `tablename` with the
name of the data frame
