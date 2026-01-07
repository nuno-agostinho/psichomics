# On collapse observer

Adds a JavaScript listener to a Bootstrap collapse panel so that when it
is shown, a Shiny input value is updated with the provided label.

## Usage

``` r
onCollapseOpen(id)
```

## Arguments

- id:

  The ID of the collapse panel.

## Value

A `tags$script` object containing the JavaScript listener.
