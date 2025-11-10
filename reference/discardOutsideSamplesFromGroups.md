# Discard grouped samples if not within a sample vector

Discard grouped samples if not within a sample vector

## Usage

``` r
discardOutsideSamplesFromGroups(groups, samples, clean = FALSE)
```

## Arguments

- groups:

  Named list of samples

- samples:

  Character: vector with all available samples

- clean:

  Boolean: clean results?

## Value

Groups without samples not found in `samples`
