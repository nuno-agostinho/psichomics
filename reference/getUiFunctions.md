# Matches user interface (UI) functions from a given loader

Matches user interface (UI) functions from a given loader

## Usage

``` r
getUiFunctions(ns, loader, ..., priority = NULL)
```

## Arguments

- ns:

  Shiny function to create IDs within a namespace

- loader:

  Character: loader to run the functions

- ...:

  Extra arguments to pass to the user interface (UI) functions

- priority:

  Character: name of functions to prioritise by the given order; for
  instance, `c("data", "analyses")` would load `data`, then `analyses`
  and finally the remaining functions

## Value

List of functions related to the given loader
