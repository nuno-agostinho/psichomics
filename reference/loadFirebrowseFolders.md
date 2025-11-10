# Load FireBrowse folders

Loads the files present in each folder as a data.frame.

## Usage

``` r
loadFirebrowseFolders(folder, exclude = "")
```

## Arguments

- folder:

  Character: folder(s) in which to look for FireBrowse files

- exclude:

  Character: files to exclude from the loading

## Value

List with loaded data.frames

## Note

For faster execution, this function uses the `readr` library. This
function ignores subfolders of the given folder (which means that files
inside subfolders are NOT loaded).
