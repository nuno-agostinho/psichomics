# Load local files

Load local files

## Usage

``` r
loadLocalFiles(
  folder,
  ignore = c(".aux.", ".mage-tab."),
  name = "Data",
  verbose = FALSE
)
```

## Arguments

- folder:

  Character: path to folder or ZIP archive

- ignore:

  Character: skip folders and filenames that match the expression

- name:

  Character: name

- verbose:

  Boolean: print steps?

## Value

List of data frames from valid files

## See also

Other functions to load data:
[`loadGtexData()`](https://nuno-agostinho.github.io/psichomics/reference/loadGtexData.md),
[`loadSRAproject()`](https://nuno-agostinho.github.io/psichomics/reference/loadSRAproject.md),
[`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
folder <- "~/Downloads/ACC 2016"
data <- loadLocalFiles(folder)

ignore <- c(".aux.", ".mage-tab.", "junction quantification")
loadLocalFiles(folder, ignore)
} # }
```
