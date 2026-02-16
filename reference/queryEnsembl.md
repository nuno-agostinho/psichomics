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
#>                            logic_name              id strand external_name
#> 1 ensembl_havana_gene_homo_sapiens_37 ENSG00000157764     -1          BRAF
#> 2               ncrna_homo_sapiens_37 ENSG00000271932      1            U6
#>   assembly_name         gene_id seq_region_name version feature_type
#> 1        GRCh37 ENSG00000157764               7       8         gene
#> 2        GRCh37 ENSG00000271932               7       1         gene
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   canonical_transcript         source       end        biotype     start
#> 1    ENST00000288602.6 ensembl_havana 140624564 protein_coding 140419127
#> 2    ENST00000605989.1        ensembl 140583978          snRNA 140583872

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $source
#> [1] "ensembl_havana"
#> 
#> $start
#> [1] 32889611
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $Transcript
#>   length         source Translation.length Translation.start
#> 1  10930 ensembl_havana               3418          32890598
#> 2   2011         havana                481          32899266
#> 3    495         havana                 64          32945108
#> 4    842         havana                186          32953977
#> 5    523         havana                 NA                NA
#> 6  10984        ensembl               3418          32890598
#>   Translation.species Translation.Parent  Translation.id Translation.end
#> 1               human    ENST00000380152 ENSP00000369497        32972907
#> 2               human    ENST00000530893 ENSP00000435699        32907428
#> 3               human    ENST00000528762 ENSP00000433168        32950807
#> 4               human    ENST00000470094 ENSP00000434898        32970229
#> 5                <NA>               <NA>            <NA>              NA
#> 6               human    ENST00000544455 ENSP00000439902        32972907
#>   Translation.version Translation.object_type Translation.db_type display_name
#> 1                   3             Translation                core    BRCA2-001
#> 2                   2             Translation                core    BRCA2-003
#> 3                   1             Translation                core    BRCA2-005
#> 4                   1             Translation                core    BRCA2-002
#> 5                  NA                    <NA>                <NA>    BRCA2-006
#> 6                   1             Translation                core    BRCA2-201
#>      start              id gencode_primary object_type                 biotype
#> 1 32889611 ENST00000380152               0  Transcript          protein_coding
#> 2 32889642 ENST00000530893               0  Transcript          protein_coding
#> 3 32945108 ENST00000528762               0  Transcript nonsense_mediated_decay
#> 4 32953977 ENST00000470094               0  Transcript nonsense_mediated_decay
#> 5 32970946 ENST00000533776               0  Transcript         retained_intron
#> 6 32889617 ENST00000544455               0  Transcript          protein_coding
#>                  logic_name version assembly_name strand species         Exon
#> 1 ensembl_havana_transcript       3        GRCh37      1   human c(328896....
#> 2    havana_homo_sapiens_37       2        GRCh37      1   human c("core"....
#> 3    havana_homo_sapiens_37       1        GRCh37      1   human c(1, 1, ....
#> 4    havana_homo_sapiens_37       1        GRCh37      1   human c(1, 1, ....
#> 5    havana_homo_sapiens_37       1        GRCh37      1   human c("Exon"....
#> 6   ensembl_homo_sapiens_37       1        GRCh37      1   human c(1, 1, ....
#>   is_canonical      end seq_region_name          Parent db_type
#> 1            0 32973347              13 ENSG00000139618    core
#> 2            0 32907428              13 ENSG00000139618    core
#> 3            0 32953632              13 ENSG00000139618    core
#> 4            0 32972409              13 ENSG00000139618    core
#> 5            0 32972585              13 ENSG00000139618    core
#> 6            1 32973805              13 ENSG00000139618    core
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $version
#> [1] 10
#> 
#> $assembly_name
#> [1] "GRCh37"
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
#> $strand
#> [1] 1
#> 
#> $species
#> [1] "human"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $end
#> [1] 32973805
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $db_type
#> [1] "core"
#> 
```
