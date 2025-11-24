# Query the Ensembl REST API

Query the Ensembl REST API

## Usage

``` r
queryEnsembl(path, query, grch37 = TRUE)
```

## Arguments

- path:

  Character: API path

- query:

  Character: API query

- grch37:

  Boolean: query the Ensembl GRCh37 API? if `FALSE`, query the most
  recent API

## Value

Parsed response or `NULL` if no response

## Examples

``` r
path  <- "overlap/region/human/7:140424943-140624564"
query <- list(feature = "gene")
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#>                            logic_name canonical_transcript       end
#> 1 ensembl_havana_gene_homo_sapiens_37    ENST00000288602.6 140624564
#> 2               ncrna_homo_sapiens_37    ENST00000605989.1 140583978
#>                id version assembly_name     start         source feature_type
#> 1 ENSG00000157764       8        GRCh37 140419127 ensembl_havana         gene
#> 2 ENSG00000271932       1        GRCh37 140583872        ensembl         gene
#>   strand external_name         gene_id seq_region_name
#> 1     -1          BRAF ENSG00000157764               7
#> 2      1            U6 ENSG00000271932               7
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>          biotype
#> 1 protein_coding
#> 2          snRNA

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $display_name
#> [1] "BRCA2"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $strand
#> [1] 1
#> 
#> $species
#> [1] "human"
#> 
#> $version
#> [1] 10
#> 
#> $end
#> [1] 32973805
#> 
#> $Transcript
#>           Exon gencode_primary                logic_name object_type length
#> 1 c("ENSE0....               0 ensembl_havana_transcript  Transcript  10930
#> 2 c(328896....               0    havana_homo_sapiens_37  Transcript   2011
#> 3 c("GRCh3....               0    havana_homo_sapiens_37  Transcript    495
#> 4 c(1, 1, ....               0    havana_homo_sapiens_37  Transcript    842
#> 5 c(329709....               0    havana_homo_sapiens_37  Transcript    523
#> 6 c("Exon"....               0   ensembl_homo_sapiens_37  Transcript  10984
#>                   biotype db_type      end    start strand species version
#> 1          protein_coding    core 32973347 32889611      1   human       3
#> 2          protein_coding    core 32907428 32889642      1   human       2
#> 3 nonsense_mediated_decay    core 32953632 32945108      1   human       1
#> 4 nonsense_mediated_decay    core 32972409 32953977      1   human       1
#> 5         retained_intron    core 32972585 32970946      1   human       1
#> 6          protein_coding    core 32973805 32889617      1   human       1
#>   is_canonical         source          Parent seq_region_name
#> 1            0 ensembl_havana ENSG00000139618              13
#> 2            0         havana ENSG00000139618              13
#> 3            0         havana ENSG00000139618              13
#> 4            0         havana ENSG00000139618              13
#> 5            0         havana ENSG00000139618              13
#> 6            1        ensembl ENSG00000139618              13
#>   Translation.species Translation.version Translation.object_type
#> 1               human                   3             Translation
#> 2               human                   2             Translation
#> 3               human                   1             Translation
#> 4               human                   1             Translation
#> 5                <NA>                  NA                    <NA>
#> 6               human                   1             Translation
#>   Translation.Parent Translation.length  Translation.id Translation.db_type
#> 1    ENST00000380152               3418 ENSP00000369497                core
#> 2    ENST00000530893                481 ENSP00000435699                core
#> 3    ENST00000528762                 64 ENSP00000433168                core
#> 4    ENST00000470094                186 ENSP00000434898                core
#> 5               <NA>                 NA            <NA>                <NA>
#> 6    ENST00000544455               3418 ENSP00000439902                core
#>   Translation.end Translation.start              id display_name assembly_name
#> 1        32972907          32890598 ENST00000380152    BRCA2-001        GRCh37
#> 2        32907428          32899266 ENST00000530893    BRCA2-003        GRCh37
#> 3        32950807          32945108 ENST00000528762    BRCA2-005        GRCh37
#> 4        32970229          32953977 ENST00000470094    BRCA2-002        GRCh37
#> 5              NA                NA ENST00000533776    BRCA2-006        GRCh37
#> 6        32972907          32890598 ENST00000544455    BRCA2-201        GRCh37
#> 
#> $db_type
#> [1] "core"
#> 
#> $start
#> [1] 32889611
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
```
