# Download files to a given directory

Download files to a given directory

## Usage

``` r
downloadFiles(url, folder, download = download.file, ...)
```

## Arguments

- url:

  Character: download links

- folder:

  Character: directory to store the downloaded archives

- download:

  Function to use to download files

- ...:

  Extra parameters passed to the download function

## Value

Invisible TRUE if every file was successfully downloaded

## Examples

``` r
if (FALSE) { # \dontrun{
url <- paste0("https://unsplash.it/400/300/?image=", 570:572)
psichomics:::downloadFiles(url, "~/Pictures")

# Download without printing to console
psichomics:::downloadFiles(url, "~/Pictures", quiet = TRUE)
} # }
```
