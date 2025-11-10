# Enable history navigation

Navigate app according to the location given by the navigation bar. Code
and logic adapted from
<https://github.com/daattali/advanced-shiny/blob/master/navigate-history>

## Usage

``` r
browserHistory(navId, input, session)
```

## Arguments

- navId:

  Character: identifier of the navigation bar

- input:

  Input object

- session:

  Session object

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
