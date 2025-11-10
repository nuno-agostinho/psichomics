# Matches server functions from a given loader

Matches server functions from a given loader

## Usage

``` r
getServerFunctions(loader, ..., priority = NULL)
```

## Arguments

- loader:

  Character: loader to run the functions

- ...:

  Extra arguments to pass to server functions

- priority:

  Character: name of functions to prioritise by the given order; for
  instance, `c("data", "analyses")` would load `data`, then `analyses`
  and finally the remaining functions

## Value

Invisible TRUE
