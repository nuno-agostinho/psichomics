# Download and load SRA projects via [recount2](https://jhubiostatistics.shinyapps.io/recount/)

Download and load SRA projects via
[recount2](https://jhubiostatistics.shinyapps.io/recount/)

## Usage

``` r
loadSRAproject(project, outdir = getDownloadsFolder())
```

## Arguments

- project:

  Character: SRA project identifiers (check
  [`recount_abstract`](https://rdrr.io/pkg/recount/man/recount_abstract.html))

- outdir:

  Character: directory to store the downloaded files

## Value

List with loaded projects

## See also

Other functions associated with SRA data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md)

Other functions to load data:
[`loadGtexData()`](https://nuno-agostinho.github.io/psichomics/reference/loadGtexData.md),
[`loadLocalFiles()`](https://nuno-agostinho.github.io/psichomics/reference/loadLocalFiles.md),
[`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
View(recount::recount_abstract)
sra <- loadSRAproject("SRP053101")
names(sra)
names(sra[[1]])
} # }
```
