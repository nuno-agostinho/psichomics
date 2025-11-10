# Missing information modal template

Missing information modal template

## Usage

``` r
loadRequiredData(modal = NULL)

missingDataModal(session, dataType, buttonId)

missingDataGuide(dataType)
```

## Arguments

- modal:

  Character: modal identifier

- session:

  Shiny session

- dataType:

  Character: type of data missing

- buttonId:

  Character: identifier of button to take user to load missing data

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## Examples

``` r
if (FALSE) { # \dontrun{
if (shiny::isRunning()) {
    session <- session$ns
    buttonInput <- "takeMeThere"
    buttonId <- ns(buttonInput)
    dataType <- "Inclusion levels"
    missingDataModal(session, buttonId, dataType)
    observeEvent(input[[buttonInput]], missingDataGuide(dataType))
}
} # }
```
