# Render a specific data tab (including data table and related interface)

Render a specific data tab (including data table and related interface)

## Usage

``` r
createDataTab(index, data, name, session, input, output)
```

## Arguments

- index:

  Integer: index of the data to load

- data:

  Data frame: data with everything to load

- name:

  Character: name of the dataset

- session:

  Shiny session

- input:

  Shiny session input

- output:

  Shiny session output

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
