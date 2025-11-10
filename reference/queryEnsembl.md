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
#>   seq_region_name         gene_id external_name assembly_name strand
#> 1               7 ENSG00000157764          BRAF        GRCh37     -1
#> 2               7 ENSG00000271932            U6        GRCh37      1
#>                            logic_name              id     start         source
#> 1 ensembl_havana_gene_homo_sapiens_37 ENSG00000157764 140419127 ensembl_havana
#> 2               ncrna_homo_sapiens_37 ENSG00000271932 140583872        ensembl
#>          biotype       end canonical_transcript version feature_type
#> 1 protein_coding 140624564    ENST00000288602.6       8         gene
#> 2          snRNA 140583978    ENST00000605989.1       1         gene
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
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
#> $db_type
#> [1] "core"
#> 
#> $Transcript
#>   version strand species is_canonical seq_region_name Translation.version
#> 1       3      1   human            0              13                   3
#> 2       2      1   human            0              13                   2
#> 3       1      1   human            0              13                   1
#> 4       1      1   human            0              13                   1
#> 5       1      1   human            0              13                  NA
#> 6       1      1   human            1              13                   1
#>   Translation.species Translation.start Translation.end Translation.db_type
#> 1               human          32890598        32972907                core
#> 2               human          32899266        32907428                core
#> 3               human          32945108        32950807                core
#> 4               human          32953977        32970229                core
#> 5                <NA>                NA              NA                <NA>
#> 6               human          32890598        32972907                core
#>    Translation.id Translation.length Translation.Parent Translation.object_type
#> 1 ENSP00000369497               3418    ENST00000380152             Translation
#> 2 ENSP00000435699                481    ENST00000530893             Translation
#> 3 ENSP00000433168                 64    ENST00000528762             Translation
#> 4 ENSP00000434898                186    ENST00000470094             Translation
#> 5            <NA>                 NA               <NA>                    <NA>
#> 6 ENSP00000439902               3418    ENST00000544455             Translation
#>                id          Parent         source assembly_name display_name
#> 1 ENST00000380152 ENSG00000139618 ensembl_havana        GRCh37    BRCA2-001
#> 2 ENST00000530893 ENSG00000139618         havana        GRCh37    BRCA2-003
#> 3 ENST00000528762 ENSG00000139618         havana        GRCh37    BRCA2-005
#> 4 ENST00000470094 ENSG00000139618         havana        GRCh37    BRCA2-002
#> 5 ENST00000533776 ENSG00000139618         havana        GRCh37    BRCA2-006
#> 6 ENST00000544455 ENSG00000139618        ensembl        GRCh37    BRCA2-201
#>           Exon gencode_primary                logic_name length
#> 1 c("human....               0 ensembl_havana_transcript  10930
#> 2 c(1, 1, ....               0    havana_homo_sapiens_37   2011
#> 3 c(1, 1, ....               0    havana_homo_sapiens_37    495
#> 4 c(1, 1, ....               0    havana_homo_sapiens_37    842
#> 5 c(329711....               0    havana_homo_sapiens_37    523
#> 6 c(1, 1, ....               0   ensembl_homo_sapiens_37  10984
#>                   biotype object_type    start db_type      end
#> 1          protein_coding  Transcript 32889611    core 32973347
#> 2          protein_coding  Transcript 32889642    core 32907428
#> 3 nonsense_mediated_decay  Transcript 32945108    core 32953632
#> 4 nonsense_mediated_decay  Transcript 32953977    core 32972409
#> 5         retained_intron  Transcript 32970946    core 32972585
#> 6          protein_coding  Transcript 32889617    core 32973805
#> 
#> $end
#> [1] 32973805
#> 
#> $start
#> [1] 32889611
#> 
#> $species
#> [1] "human"
#> 
#> $strand
#> [1] 1
#> 
#> $version
#> [1] 10
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
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
#> $display_name
#> [1] "BRCA2"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
```
