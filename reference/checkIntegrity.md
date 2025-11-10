# Compute the 32-byte `MD5` hashes of one or more files and check with given `md5` file

Compute the 32-byte `MD5` hashes of one or more files and check with
given `md5` file

## Usage

``` r
checkIntegrity(filesToCheck, md5file)
```

## Arguments

- filesToCheck:

  Character: files to calculate and match `MD5` hashes

- md5file:

  Character: file containing correct `MD5` hashes

## Value

Logical vector showing TRUE for files with matching `md5sums` and
`FALSE` for files with non-matching `md5sums`
