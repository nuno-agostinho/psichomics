# Convert gene identifiers

Convert gene identifiers

## Usage

``` r
convertGeneIdentifiers(
  annotation,
  genes,
  key = "ENSEMBL",
  target = "SYMBOL",
  ignoreDuplicatedTargets = TRUE
)
```

## Arguments

- annotation:

  `OrgDb` with genome wide annotation for an organism or `character`
  with species name to query `OrgDb`, e.g. `"Homo sapiens"`

- genes:

  Character: genes to be converted

- key:

  Character: type of identifier used, e.g. `ENSEMBL`; read
  [`?AnnotationDbi::columns`](https://rdrr.io/pkg/AnnotationDbi/man/AnnotationDb-class.html)

- target:

  Character: type of identifier to convert to; read
  [`?AnnotationDbi::columns`](https://rdrr.io/pkg/AnnotationDbi/man/AnnotationDb-class.html)

- ignoreDuplicatedTargets:

  Boolean: if `TRUE`, identifiers that share targets with other
  identifiers will not be converted

## Value

Character vector of the respective targets of gene identifiers. The
previous identifiers remain other identifiers have the same target (in
case `ignoreDuplicatedTargets = TRUE`) or if no target was found.

## See also

Other functions for gene expression pre-processing:
[`filterGeneExpr()`](https://nuno-agostinho.github.io/psichomics/reference/filterGeneExpr.md),
[`normaliseGeneExpression()`](https://nuno-agostinho.github.io/psichomics/reference/normaliseGeneExpression.md),
[`plotGeneExprPerSample()`](https://nuno-agostinho.github.io/psichomics/reference/plotGeneExprPerSample.md),
[`plotLibrarySize()`](https://nuno-agostinho.github.io/psichomics/reference/plotLibrarySize.md),
[`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md)

## Examples

``` r
# Use species name to automatically look for a OrgDb database
sp <- "Homo sapiens"
genes <- c("ENSG00000012048", "ENSG00000083093", "ENSG00000141510",
           "ENSG00000051180")
convertGeneIdentifiers(sp, genes)
#> loading from cache
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following object is masked from ‘package:psichomics’:
#> 
#>     plotPCA
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> ENSG00000012048 ENSG00000083093 ENSG00000141510 ENSG00000051180 
#>         "BRCA1"         "PALB2"          "TP53"         "RAD51" 
convertGeneIdentifiers(sp, genes, key="ENSEMBL", target="UNIPROT")
#> loading from cache
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000012048 
#>       "A0A386INP3/A0A386IPK6/A0A2R8Y7V5/G8I0D8/E9PFZ0/O15129/P38398/Q1RMC1/Q3LRJ0/Q3LRJ6/Q6IN79/Q7KYU9/H0Y850/E7EQW4/A0A9Y1VUM8/E7ENB7/E9PH68/A0A9Y1QQD3/A0A0U1RRA9/H0Y8D8/A0A9Y1QQJ6/A0A386IN41/A0A9Y1QQK3/A0A9Y1QQK7/A0A386IN52/A0A9Y1QPY6/A0A9Y1VR53/C9IZW4/A0A9Y1QQK5/A0A9Y1QQL7/A0A9Y1VVE2/E7EUM2/A0A9Y1VVF5/A0A9Y1VVD0/H0Y8B8/A0A9Y1QPR4/A0A9Y1QQF1/A0A2R8Y587/A0A9Y1QPT7/B4DES0/A0A9Y1QQ47/A0A9Y1QPQ7/A0A9Y1QQ02/A0A9Y1VVF6/A0A9Y1QQ22" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000083093 
#>                                                                                                                                                                                                                                                                                                  "A0AA52I2C1/A0A8V8TKZ4/A0AAQ5BGU0/A0AAQ5BGW2/I3L1Z5/B4DR89/H3BN63/I3L3R6/A0A8V8TMC9/A0A8V8TMK8/A0A8V8TLC8/A6NIE1/Q86YC2/Q8N7Y6/Q8ND31/Q9H6W1" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000141510 
#> "A0A386NC20/A0A386NC22/A0A386NC55/A0A386NC62/A0A386NCA6/A0A386NCB1/A0A386NCU2/A0A386NCW8/A0A386NCX4/A0A386NDA8/A0A386NDB3/A0A386NFX3/A0A386NFY2/A0A386NG45/K7PPA8/K7PPU4/P04637/Q15086/Q15087/Q15088/Q16535/Q16807/Q16808/Q16809/Q16810/Q16811/Q16848/Q2XN98/Q3LRW1/Q3LRW2/Q3LRW3/Q3LRW4/Q3LRW5/Q86UG1/Q8J016/Q99659/Q9BTM4/Q9HAQ8/Q9NP68/Q9NPJ2/Q9NZD0/Q9UBI2/Q9UQ61/A0A223PQI5/J3KP33/E9PFT5/E7ESS1/H2EHT1/A0A087X1Q1/A0A087WXZ1/A0A087WT22" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000051180 
#>                                                                                                                                                                                                                                                                                                                                                                                             "Q5U0A5/Q06609/B0FXP0/B2R8T6/Q6FHX9/Q6ZNA8/Q9BV60" 

# Alternatively, set the annotation database directly
ah <- AnnotationHub::AnnotationHub()
sp <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]
#> loading from cache
columns(sp) # these attributes can be used to change the attributes
#>  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
#>  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
#> [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
#> [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
#> [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
#> [26] "UNIPROT"     

convertGeneIdentifiers(sp, genes)
#> ENSG00000012048 ENSG00000083093 ENSG00000141510 ENSG00000051180 
#>         "BRCA1"         "PALB2"          "TP53"         "RAD51" 
convertGeneIdentifiers(sp, genes, key="ENSEMBL", target="UNIPROT")
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000012048 
#>       "A0A386INP3/A0A386IPK6/A0A2R8Y7V5/G8I0D8/E9PFZ0/O15129/P38398/Q1RMC1/Q3LRJ0/Q3LRJ6/Q6IN79/Q7KYU9/H0Y850/E7EQW4/A0A9Y1VUM8/E7ENB7/E9PH68/A0A9Y1QQD3/A0A0U1RRA9/H0Y8D8/A0A9Y1QQJ6/A0A386IN41/A0A9Y1QQK3/A0A9Y1QQK7/A0A386IN52/A0A9Y1QPY6/A0A9Y1VR53/C9IZW4/A0A9Y1QQK5/A0A9Y1QQL7/A0A9Y1VVE2/E7EUM2/A0A9Y1VVF5/A0A9Y1VVD0/H0Y8B8/A0A9Y1QPR4/A0A9Y1QQF1/A0A2R8Y587/A0A9Y1QPT7/B4DES0/A0A9Y1QQ47/A0A9Y1QPQ7/A0A9Y1QQ02/A0A9Y1VVF6/A0A9Y1QQ22" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000083093 
#>                                                                                                                                                                                                                                                                                                  "A0AA52I2C1/A0A8V8TKZ4/A0AAQ5BGU0/A0AAQ5BGW2/I3L1Z5/B4DR89/H3BN63/I3L3R6/A0A8V8TMC9/A0A8V8TMK8/A0A8V8TLC8/A6NIE1/Q86YC2/Q8N7Y6/Q8ND31/Q9H6W1" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000141510 
#> "A0A386NC20/A0A386NC22/A0A386NC55/A0A386NC62/A0A386NCA6/A0A386NCB1/A0A386NCU2/A0A386NCW8/A0A386NCX4/A0A386NDA8/A0A386NDB3/A0A386NFX3/A0A386NFY2/A0A386NG45/K7PPA8/K7PPU4/P04637/Q15086/Q15087/Q15088/Q16535/Q16807/Q16808/Q16809/Q16810/Q16811/Q16848/Q2XN98/Q3LRW1/Q3LRW2/Q3LRW3/Q3LRW4/Q3LRW5/Q86UG1/Q8J016/Q99659/Q9BTM4/Q9HAQ8/Q9NP68/Q9NPJ2/Q9NZD0/Q9UBI2/Q9UQ61/A0A223PQI5/J3KP33/E9PFT5/E7ESS1/H2EHT1/A0A087X1Q1/A0A087WXZ1/A0A087WT22" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                ENSG00000051180 
#>                                                                                                                                                                                                                                                                                                                                                                                             "Q5U0A5/Q06609/B0FXP0/B2R8T6/Q6FHX9/Q6ZNA8/Q9BV60" 
```
