# Return the interface of relevant PubMed articles for a given gene

Return the interface of relevant PubMed articles for a given gene

## Usage

``` r
pubmedUI(ns, gene, ...)
```

## Arguments

- ns:

  Namespace function

- gene:

  Character: gene

- ...:

  Arguments passed on to
  [`queryPubMed`](https://nuno-agostinho.github.io/psichomics/reference/queryPubMed.md)

  `top`

  :   Numeric: number of articles to retrieve

  `field`

  :   Character: field of interest where to look for terms (`abstract`
      by default)

  `sort`

  :   Character: sort by a given parameter (`relevance` by default)

## Value

HTML interface of relevant PubMed articles
