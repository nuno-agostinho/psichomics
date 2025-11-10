# Preserve attributes of `sticky` objects when extracting or transposing object

Most attributes - with the exception of `names`, `dim`, `dimnames`,
`class` and `row.names` - are preserved in simple transformations of
objects from class `sticky`

## Usage

``` r
# S3 method for class 'sticky'
t(x)

# S3 method for class 'sticky'
x[i, j, ...]
```

## Arguments

- x:

  Object

- i, j, ...:

  Numeric or character: indices of elements to extract

## Value

Transformed object with most attributes preserved
