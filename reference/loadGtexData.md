# Download and load GTEx data

Download and load GTEx data

## Usage

``` r
loadGtexData(
  folder = getDownloadsFolder(),
  data = getGtexDataTypes(),
  tissue = NULL,
  release = getGtexReleases()[[1]],
  progress = TRUE
)
```

## Arguments

- folder:

  Character: folder containing data

- data:

  Character: data types to load (see `getGtexDataTypes`)

- tissue:

  Character: tissues to load (if `NULL`, load all); tissue selection may
  speed up data loading

- release:

  Numeric: GTEx data release to load

- progress:

  Boolean: display progress?

## Value

List with loaded data

## See also

Other functions associated with GTEx data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`getGtexDataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexDataTypes.md),
[`getGtexTissues()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexTissues.md)

Other functions to load data:
[`loadLocalFiles()`](https://nuno-agostinho.github.io/psichomics/reference/loadLocalFiles.md),
[`loadSRAproject()`](https://nuno-agostinho.github.io/psichomics/reference/loadSRAproject.md),
[`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Download and load all available GTEx data
data <- loadGtexData()

# Download and load only junction quantification and sample info from GTEx
getGtexDataTypes()
data <- loadGtexData(data=c("sampleInfo", "junctionQuant"))

# Download and load only data for specific tissues
getGtexTissues()
data <- loadGtexData(tissue=c("Stomach", "Small Intestine"))

# Download and load data from a specific GTEx data release
data <- loadGtexData(tissue=c("Stomach", "Small Intestine"), release=7)
} # }
```
