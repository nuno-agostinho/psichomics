# Create word break opportunities (for HTML) using given characters

Create word break opportunities (for HTML) using given characters

## Usage

``` r
prepareWordBreak(
  str,
  pattern = c(".", "-", "\\", "/", "_", ",", " ", "+", "="),
  html = TRUE
)
```

## Arguments

- str:

  Character: text

- pattern:

  Character: pattern(s) of interest to be used as word break
  opportunities

- html:

  Boolean: convert to HTML?

## Value

String containing HTML elements
