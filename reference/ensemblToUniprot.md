# Convert from Ensembl to UniProt identifier

Convert from Ensembl to UniProt identifier

## Usage

``` r
ensemblToUniprot(protein)
```

## Arguments

- protein:

  Character: Ensembl identifier

## Value

UniProt protein identifier

## See also

Other functions to retrieve external information:
[`plotProtein()`](https://nuno-agostinho.github.io/psichomics/reference/plotProtein.md),
[`plotTranscripts()`](https://nuno-agostinho.github.io/psichomics/reference/plotTranscripts.md),
[`queryEnsemblByGene()`](https://nuno-agostinho.github.io/psichomics/reference/queryEnsemblByGene.md)

## Examples

``` r
gene <- "ENSG00000173262"
ensemblToUniprot(gene)
#> SLC2A14 (UniProtKB Gene Name) 
#>                     "SLC2A14" 

protein <- "ENSP00000445929"
ensemblToUniprot(protein)
#> B7ZAC3_HUMAN (UniProtKB/TrEMBL) 
#>                        "B7ZAC3" 
```
