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
#>           source feature_type         gene_id
#> 1 ensembl_havana         gene ENSG00000157764
#> 2        ensembl         gene ENSG00000271932
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>                id                          logic_name        biotype
#> 1 ENSG00000157764 ensembl_havana_gene_homo_sapiens_37 protein_coding
#> 2 ENSG00000271932               ncrna_homo_sapiens_37          snRNA
#>   external_name strand       end canonical_transcript version seq_region_name
#> 1          BRAF     -1 140624564    ENST00000288602.6       8               7
#> 2            U6      1 140583978    ENST00000605989.1       1               7
#>   assembly_name     start
#> 1        GRCh37 140419127
#> 2        GRCh37 140583872

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $display_name
#> [1] "BRCA2"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $start
#> [1] 32889611
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $version
#> [1] 10
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $Transcript
#>                   biotype length    start version          Parent object_type
#> 1          protein_coding  10930 32889611       3 ENSG00000139618  Transcript
#> 2          protein_coding   2011 32889642       2 ENSG00000139618  Transcript
#> 3 nonsense_mediated_decay    495 32945108       1 ENSG00000139618  Transcript
#> 4 nonsense_mediated_decay    842 32953977       1 ENSG00000139618  Transcript
#> 5         retained_intron    523 32970946       1 ENSG00000139618  Transcript
#> 6          protein_coding  10984 32889617       1 ENSG00000139618  Transcript
#>   display_name  Translation.id Translation.db_type Translation.object_type
#> 1    BRCA2-001 ENSP00000369497                core             Translation
#> 2    BRCA2-003 ENSP00000435699                core             Translation
#> 3    BRCA2-005 ENSP00000433168                core             Translation
#> 4    BRCA2-002 ENSP00000434898                core             Translation
#> 5    BRCA2-006            <NA>                <NA>                    <NA>
#> 6    BRCA2-201 ENSP00000439902                core             Translation
#>   Translation.species Translation.Parent Translation.version Translation.start
#> 1               human    ENST00000380152                   3          32890598
#> 2               human    ENST00000530893                   2          32899266
#> 3               human    ENST00000528762                   1          32945108
#> 4               human    ENST00000470094                   1          32953977
#> 5                <NA>               <NA>                  NA                NA
#> 6               human    ENST00000544455                   1          32890598
#>   Translation.length Translation.end         Exon is_canonical      end
#> 1               3418        32972907 c(1, 1, ....            0 32973347
#> 2                481        32907428 c(1, 1, ....            0 32907428
#> 3                 64        32950807 c(1, 1, ....            0 32953632
#> 4                186        32970229 c(329540....            0 32972409
#> 5                 NA              NA c("GRCh3....            0 32972585
#> 6               3418        32972907 c("GRCh3....            1 32973805
#>           source assembly_name species seq_region_name              id
#> 1 ensembl_havana        GRCh37   human              13 ENST00000380152
#> 2         havana        GRCh37   human              13 ENST00000530893
#> 3         havana        GRCh37   human              13 ENST00000528762
#> 4         havana        GRCh37   human              13 ENST00000470094
#> 5         havana        GRCh37   human              13 ENST00000533776
#> 6        ensembl        GRCh37   human              13 ENST00000544455
#>   gencode_primary db_type                logic_name strand
#> 1               0    core ensembl_havana_transcript      1
#> 2               0    core    havana_homo_sapiens_37      1
#> 3               0    core    havana_homo_sapiens_37      1
#> 4               0    core    havana_homo_sapiens_37      1
#> 5               0    core    havana_homo_sapiens_37      1
#> 6               0    core   ensembl_homo_sapiens_37      1
#> 
#> $species
#> [1] "human"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $strand
#> [1] 1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $db_type
#> [1] "core"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $end
#> [1] 32973805
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
```
