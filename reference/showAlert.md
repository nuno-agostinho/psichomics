# Show or remove an alert

Show or remove an alert

## Usage

``` r
showAlert(
  session,
  ...,
  title,
  style = NULL,
  dismissible = TRUE,
  alertId = "alert",
  iconName = NULL,
  caller = NULL
)

successAlert(
  session,
  ...,
  title = NULL,
  dismissible = TRUE,
  alertId = "success",
  caller = NULL
)

errorAlert(
  session,
  ...,
  title = NULL,
  dismissible = TRUE,
  alertId = "alert",
  caller = NULL
)

warningAlert(
  session,
  ...,
  title = NULL,
  dismissible = TRUE,
  alertId = "alert",
  caller = NULL
)

removeAlert(output, alertId = "alert")
```

## Arguments

- session:

  Shiny session

- ...:

  Arguments to render as elements of alert

- title:

  Character: title

- style:

  Character: style (`error`, `warning` or `NULL`)

- dismissible:

  Boolean: is the alert dismissible?

- alertId:

  Character: identifier

- iconName:

  Character: icon name

- caller:

  Character: caller module identifier

- output:

  Shiny output

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## See also

[`showModal()`](https://rdrr.io/pkg/shiny/man/showModal.html)
