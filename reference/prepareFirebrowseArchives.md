# Prepares FireBrowse archives in a given directory

Checks FireBrowse archives' integrity using the MD5 files, extracts the
content of the archives, moves the content to newly-created folders and
removes the original downloaded archives.

## Usage

``` r
prepareFirebrowseArchives(archive, md5, folder, outdir)
```

## Arguments

- archive:

  Character: path to downloaded archives

- md5:

  Character: path to MD5 files of each archive

- folder:

  Character: master directory where every archive will be extracted

- outdir:

  Character: subdirectories where to move the extracted content

## Value

Invisible TRUE if successful

## Examples

``` r
file <- paste0(
    "~/Downloads",
    "ACC/20151101/gdac.broadinstitute.org_ACC.",
    "Merge_Clinical.Level_1.2015110100.0.0.tar.gz")
md5 <- paste0(file, ".md5")
if (FALSE) { # \dontrun{
prepareFirebrowseArchives(archive = file, md5 = paste0(file, ".md5"))
} # }
```
