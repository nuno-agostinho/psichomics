# Prepare file browser dialogue and update the input's value accordingly to selected file or directory

Prepare file browser dialogue and update the input's value accordingly
to selected file or directory

## Usage

``` r
prepareFileBrowser(session, input, id, modalId = "modal", ...)
```

## Arguments

- session:

  Shiny session

- input:

  Shiny input

- id:

  Character: input identifier

- modalId:

  Character: modal window identifier

- ...:

  Arguments passed on to
  [`fileBrowser`](https://nuno-agostinho.github.io/psichomics/reference/fileBrowser.md)

  `default`

  :   Character: path to initial folder

  `caption`

  :   Character: caption on the selection dialogue

  `multiple`

  :   Boolean: allow to select multiple files?

  `directory`

  :   Boolean: allow to select directories instead of files?

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
