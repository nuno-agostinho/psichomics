# Look for event data in input

Check if event data can be found in `data` and then `event`. Event data
has to be an object of class `eventData`

## Usage

``` r
findEventData(event = NULL, data = NULL)
```

## Arguments

- event:

  Character: AS event that may contain event data in its attribute
  `eventData`

- data:

  Data frame or matrix: either event data or data containing event data
  in its attributes `rowData` or `eventData`

## Value

Event data (or `NULL` if not found)
