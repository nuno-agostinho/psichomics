# Get rows of a data frame between two row indexes

Get rows of a data frame between two row indexes

## Usage

``` r
getDataRows(i, data, firstRow, lastRow)
```

## Arguments

- i:

  Integer: current iteration

- data:

  Data.frame: contains the data of interest

- firstRow:

  Vector of integers: First row index of interest; value must be less
  than the respective last row index and less than the number of rows in
  the data frame

- lastRow:

  Vector of integers: Last row index of interest; value must be higher
  than the respective first row index and less than the number of rows
  in the data frame

## Value

Data frame subset from two row indexes (returns NA if the first row
index is NA)

## Details

For a given iteration i, returns data from `firstRow[i]` to `lastRow[i]`
