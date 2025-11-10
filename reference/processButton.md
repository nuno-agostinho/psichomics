# Style button used to initiate a process

Style button used to initiate a process

## Usage

``` r
processButton(id, label, ..., class = "btn-primary")
```

## Arguments

- id:

  Character: button identifier

- label:

  Character: label

- ...:

  Arguments passed on to
  [`shiny::actionButton`](https://rdrr.io/pkg/shiny/man/actionButton.html)

  `icon`

  :   An optional [`icon()`](https://rdrr.io/pkg/shiny/man/icon.html) to
      appear on the button.

  `width`

  :   The width of the input, e.g. `'400px'`, or `'100%'`; see
      [`validateCssUnit()`](https://rstudio.github.io/htmltools/reference/validateCssUnit.html).

- class:

  Character: class

## Value

HTML for a button
