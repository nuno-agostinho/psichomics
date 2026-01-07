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
#>           gene_id external_name        biotype
#> 1 ENSG00000157764          BRAF protein_coding
#> 2 ENSG00000271932            U6          snRNA
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   seq_region_name assembly_name version              id       end
#> 1               7        GRCh37       8 ENSG00000157764 140624564
#> 2               7        GRCh37       1 ENSG00000271932 140583978
#>   canonical_transcript                          logic_name strand
#> 1    ENST00000288602.6 ensembl_havana_gene_homo_sapiens_37     -1
#> 2    ENST00000605989.1               ncrna_homo_sapiens_37      1
#>           source feature_type     start
#> 1 ensembl_havana         gene 140419127
#> 2        ensembl         gene 140583872

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $display_name
#> [1] "BRCA2"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $version
#> [1] 10
#> 
#> $start
#> [1] 32889611
#> 
#> $biotype
#> [1] "protein_coding"
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
#> $Transcript
#>           source      end                logic_name gencode_primary
#> 1 ensembl_havana 32973347 ensembl_havana_transcript               0
#> 2         havana 32907428    havana_homo_sapiens_37               0
#> 3         havana 32953632    havana_homo_sapiens_37               0
#> 4         havana 32972409    havana_homo_sapiens_37               0
#> 5         havana 32972585    havana_homo_sapiens_37               0
#> 6        ensembl 32973805   ensembl_homo_sapiens_37               0
#>                id db_type strand assembly_name seq_region_name species version
#> 1 ENST00000380152    core      1        GRCh37              13   human       3
#> 2 ENST00000530893    core      1        GRCh37              13   human       2
#> 3 ENST00000528762    core      1        GRCh37              13   human       1
#> 4 ENST00000470094    core      1        GRCh37              13   human       1
#> 5 ENST00000533776    core      1        GRCh37              13   human       1
#> 6 ENST00000544455    core      1        GRCh37              13   human       1
#>   length                 biotype    start Translation.species
#> 1  10930          protein_coding 32889611               human
#> 2   2011          protein_coding 32889642               human
#> 3    495 nonsense_mediated_decay 32945108               human
#> 4    842 nonsense_mediated_decay 32953977               human
#> 5    523         retained_intron 32970946                <NA>
#> 6  10984          protein_coding 32889617               human
#>   Translation.Parent Translation.db_type Translation.object_type
#> 1    ENST00000380152                core             Translation
#> 2    ENST00000530893                core             Translation
#> 3    ENST00000528762                core             Translation
#> 4    ENST00000470094                core             Translation
#> 5               <NA>                <NA>                    <NA>
#> 6    ENST00000544455                core             Translation
#>    Translation.id Translation.start Translation.end Translation.length
#> 1 ENSP00000369497          32890598        32972907               3418
#> 2 ENSP00000435699          32899266        32907428                481
#> 3 ENSP00000433168          32945108        32950807                 64
#> 4 ENSP00000434898          32953977        32970229                186
#> 5            <NA>                NA              NA                 NA
#> 6 ENSP00000439902          32890598        32972907               3418
#>   Translation.version display_name object_type is_canonical         Exon
#> 1                   3    BRCA2-001  Transcript            0 c("GRCh3....
#> 2                   2    BRCA2-003  Transcript            0 c(1, 1, ....
#> 3                   1    BRCA2-005  Transcript            0 c(1, 1, ....
#> 4                   1    BRCA2-002  Transcript            0 c(329540....
#> 5                  NA    BRCA2-006  Transcript            0 c(1, 1),....
#> 6                   1    BRCA2-201  Transcript            1 c(1, 1, ....
#>            Parent
#> 1 ENSG00000139618
#> 2 ENSG00000139618
#> 3 ENSG00000139618
#> 4 ENSG00000139618
#> 5 ENSG00000139618
#> 6 ENSG00000139618
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $species
#> [1] "human"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $end
#> [1] 32973805
#> 
```
