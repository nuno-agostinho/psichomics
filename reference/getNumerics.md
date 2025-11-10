# Convert a column to numeric if possible and ignore given columns composed of lists

Convert a column to numeric if possible and ignore given columns
composed of lists

## Usage

``` r
getNumerics(table, by = NULL, toNumeric = FALSE)
```

## Arguments

- table:

  Data matrix: table

- by:

  Character: column names of interest

- toNumeric:

  Boolean: which columns to convert to numeric

## Value

Processed data matrix

## Examples

``` r
event <- read.table(text = "ABC123 + 250 300 350
                            DEF456 - 900 800 700")
names(event) <- c("Event ID", "Strand", "C1.end", "A1.end", "A1.start")

# Let's change one column to character
event[ , "C1.end"] <- as.character(event[ , "C1.end"])
is.character(event[ , "C1.end"])
#> [1] TRUE

event <- psichomics:::getNumerics(event, by = c("Strand", "C1.end", "A1.end",
                                  "A1.start"),
                                  toNumeric = c(FALSE, TRUE, TRUE, TRUE))
# Let's check if the same column is now integer
is.numeric(event[ , "C1.end"])
#> [1] TRUE
```
