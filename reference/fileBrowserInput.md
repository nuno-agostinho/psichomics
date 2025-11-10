# File browser input

Input to interactively select a file or directory on the server

## Usage

``` r
fileBrowserInfoInput(id, label, infoContent = NULL, clearable = FALSE)

fileBrowserInput(
  id,
  label,
  value = NULL,
  placeholder = NULL,
  info = FALSE,
  infoFUN = NULL,
  infoPlacement = "right",
  infoTitle = "",
  infoContent = "",
  clearable = FALSE
)
```

## Source

<https://github.com/wleepang/shiny-directory-input>

## Arguments

- id:

  Character: input identifier

- label:

  Character: input label (if `NULL`, no labels are displayed)

- infoContent:

  Character: text to show as content of information

- clearable:

  Boolean: allow to clear selected file or directory?

- value:

  Character: initial value (paths are expanded via
  [`path.expand()`](https://rdrr.io/r/base/path.expand.html))

- placeholder:

  Character: placeholder when no file or folder is selected

- info:

  Boolean: add information icon for tooltips and pop-overs

- infoFUN:

  Function to use to provide information (e.g.
  [`shinyBS::bsTooltip`](https://rdrr.io/pkg/shinyBS/man/bsTooltip.html)
  and
  [`shinyBS::bsPopover`](https://rdrr.io/pkg/shinyBS/man/bsPopover.html))

- infoPlacement:

  Character: placement of the information (top, bottom, right or left)

- infoTitle:

  Character: text to show as title of information

## Value

HTML elements for a file browser input

## Details

To show the dialog for file input, the
[`prepareFileBrowser()`](https://nuno-agostinho.github.io/psichomics/reference/prepareFileBrowser.md)
function needs to be included in the server logic.

This widget relies on
[`fileBrowser()`](https://nuno-agostinho.github.io/psichomics/reference/fileBrowser.md)
to present an interactive dialogue to users for selecting a directory on
the local filesystem. Therefore, this widget is intended for shiny apps
that are run locally - i.e. on the same system that files/directories
are to be accessed - and not from hosted applications (e.g. from
<https://www.shinyapps.io>).

## See also

[`updateFileBrowserInput()`](https://nuno-agostinho.github.io/psichomics/reference/updateFileBrowserInput.md)
and
[`prepareFileBrowser()`](https://nuno-agostinho.github.io/psichomics/reference/prepareFileBrowser.md)
