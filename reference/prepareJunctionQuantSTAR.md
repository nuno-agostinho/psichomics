# Prepare user-provided files to be loaded into psichomics

Prepare user-provided files to be loaded into psichomics

## Usage

``` r
prepareJunctionQuantSTAR(..., startOffset = -1, endOffset = +1)

prepareGeneQuantSTAR(
  ...,
  strandedness = c("unstranded", "stranded", "stranded (reverse)")
)
```

## Arguments

- ...:

  Character: path of (optionally named) input files (see Examples)

- startOffset:

  Numeric: value to offset start position

- endOffset:

  Numeric: value to offset end position

- strandedness:

  Character: strandedness of RNA-seq protocol; may be one of the
  following: `unstraded`, `stranded` or `stranded (reverse)`

## Value

Prepared file (if `output != NULL`) and object

## Examples

``` r
if (FALSE) { # \dontrun{
prepareJunctionQuant("Control rep1"=junctionFile1,
                     "Control rep2"=junctionFile2,
                     "KD rep1"=junctionFile3,
                     "KD rep2"=junctionFile4)
} # }
if (FALSE) { # \dontrun{
prepareGeneQuant("Control rep1"=geneCountFile1,
                 "Control rep2"=geneCountFile2,
                 "KD rep1"=geneCountFile3,
                 "KD rep2"=geneCountFile4)
} # }
```
