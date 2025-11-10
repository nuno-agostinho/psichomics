# Interactive folder selection using a native dialogue

Interactive folder selection using a native dialogue

## Usage

``` r
fileBrowser(
  default = NULL,
  caption = NULL,
  multiple = FALSE,
  directory = FALSE
)
```

## Source

<https://github.com/wleepang/shiny-directory-input>

## Arguments

- default:

  Character: path to initial folder

- caption:

  Character: caption on the selection dialogue

- multiple:

  Boolean: allow to select multiple files?

- directory:

  Boolean: allow to select directories instead of files?

## Value

A length one character vector, character NA if 'Cancel' was selected

## Details

Platform-dependent implementation:

- **Windows**: calls the
  [`utils::choose.files`](https://rdrr.io/r/utils/choose.files.html) R
  function.

- **macOS**: uses AppleScript to display a folder selection dialogue. If
  `default = NA`, folder selection falls back to the default behaviour
  of the `choose folder` AppleScript command. Otherwise, paths are
  expanded with
  [`path.expand()`](https://rdrr.io/r/base/path.expand.html).

- **Linux**: calls the `zenity` system command.
