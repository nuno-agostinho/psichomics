# Create a modal window

Create a modal window

## Usage

``` r
styleModal(
  session,
  title,
  ...,
  style = NULL,
  iconName = "exclamation-circle",
  footer = NULL,
  echo = FALSE,
  size = "medium",
  dismissButton = TRUE,
  caller = NULL
)

errorModal(session, title, ..., size = "small", footer = NULL, caller = NULL)

warningModal(session, title, ..., size = "small", footer = NULL, caller = NULL)

infoModal(session, title, ..., size = "small", footer = NULL, caller = NULL)
```

## Arguments

- session:

  Shiny session

- title:

  Character: title

- ...:

  Arguments passed on to
  [`shiny::modalDialog`](https://rdrr.io/pkg/shiny/man/modalDialog.html)

  `easyClose`

  :   If `TRUE`, the modal dialog can be dismissed by clicking outside
      the dialog box, or be pressing the Escape key. If `FALSE` (the
      default), the modal dialog can't be dismissed in those ways;
      instead it must be dismissed by clicking on a
      [`modalButton()`](https://rdrr.io/pkg/shiny/man/modalDialog.html),
      or from a call to
      [`removeModal()`](https://rdrr.io/pkg/shiny/man/showModal.html) on
      the server.

  `fade`

  :   If `FALSE`, the modal dialog will have no fade-in animation (it
      will simply appear rather than fade in to view).

- style:

  Character: style (`NULL`, `warning`, `error` or `info`)

- iconName:

  Character: icon name

- footer:

  HTML elements to use in footer

- echo:

  Boolean: print to console?

- size:

  Character: size of the modal (`small`, `medium` or `large`)

- dismissButton:

  Boolean: show dismiss button in footer?

- caller:

  Character: caller module identifier

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## See also

[`showAlert()`](https://nuno-agostinho.github.io/psichomics/reference/showAlert.md)
