# Render a data table with sparkline HTML elements

Render a data table with sparkline HTML elements

## Usage

``` r
renderDataTableSparklines(..., options = NULL)
```

## Arguments

- ...:

  Arguments passed on to
  [`shiny::renderDataTable`](https://rdrr.io/pkg/shiny/man/renderDataTable.html)

  `expr`

  :   An expression that returns a data frame or a matrix.

  `searchDelay`

  :   The delay for searching, in milliseconds (to avoid too frequent
      search requests).

  `callback`

  :   A JavaScript function to be applied to the DataTable object. This
      is useful for DataTables plug-ins, which often require the
      DataTable instance to be available.

  `quoted`

  :   If it is `TRUE`, then the
      [`quote()`](https://rdrr.io/r/base/substitute.html)ed value of
      `expr` will be used when `expr` is evaluated. If `expr` is a
      quosure and you would like to use its expression as a value for
      `expr`, then you must set `quoted` to `TRUE`.

  `outputArgs`

  :   A list of arguments to be passed through to the implicit call to
      [`dataTableOutput()`](https://rdrr.io/pkg/shiny/man/renderDataTable.html)
      when
      [`renderDataTable()`](https://rdrr.io/pkg/shiny/man/renderDataTable.html)
      is used in an interactive R Markdown document.

- options:

  List of options to pass to
  [`renderDataTable()`](https://rdrr.io/pkg/shiny/man/renderDataTable.html)

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## Details

This slightly modified version of
[`renderDataTable()`](https://rdrr.io/pkg/shiny/man/renderDataTable.html)
calls a JavaScript function to convert the sparkline HTML elements to an
interactive `highchart` object
