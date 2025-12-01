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
#>           gene_id external_name canonical_transcript seq_region_name
#> 1 ENSG00000157764          BRAF    ENST00000288602.6               7
#> 2 ENSG00000271932            U6    ENST00000605989.1               7
#>   assembly_name        biotype version
#> 1        GRCh37 protein_coding       8
#> 2        GRCh37          snRNA       1
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   feature_type     start       end         source              id
#> 1         gene 140419127 140624564 ensembl_havana ENSG00000157764
#> 2         gene 140583872 140583978        ensembl ENSG00000271932
#>                            logic_name strand
#> 1 ensembl_havana_gene_homo_sapiens_37     -1
#> 2               ncrna_homo_sapiens_37      1

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $species
#> [1] "human"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $end
#> [1] 32973805
#> 
#> $version
#> [1] 10
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $db_type
#> [1] "core"
#> 
#> $start
#> [1] 32889611
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $Transcript
#>   species seq_region_name          Parent assembly_name db_type version
#> 1   human              13 ENSG00000139618        GRCh37    core       3
#> 2   human              13 ENSG00000139618        GRCh37    core       2
#> 3   human              13 ENSG00000139618        GRCh37    core       1
#> 4   human              13 ENSG00000139618        GRCh37    core       1
#> 5   human              13 ENSG00000139618        GRCh37    core       1
#> 6   human              13 ENSG00000139618        GRCh37    core       1
#>        end Translation.object_type Translation.species Translation.Parent
#> 1 32973347             Translation               human    ENST00000380152
#> 2 32907428             Translation               human    ENST00000530893
#> 3 32953632             Translation               human    ENST00000528762
#> 4 32972409             Translation               human    ENST00000470094
#> 5 32972585                    <NA>                <NA>               <NA>
#> 6 32973805             Translation               human    ENST00000544455
#>   Translation.end  Translation.id Translation.version Translation.db_type
#> 1        32972907 ENSP00000369497                   3                core
#> 2        32907428 ENSP00000435699                   2                core
#> 3        32950807 ENSP00000433168                   1                core
#> 4        32970229 ENSP00000434898                   1                core
#> 5              NA            <NA>                  NA                <NA>
#> 6        32972907 ENSP00000439902                   1                core
#>   Translation.start Translation.length gencode_primary display_name length
#> 1          32890598               3418               0    BRCA2-001  10930
#> 2          32899266                481               0    BRCA2-003   2011
#> 3          32945108                 64               0    BRCA2-005    495
#> 4          32953977                186               0    BRCA2-002    842
#> 5                NA                 NA               0    BRCA2-006    523
#> 6          32890598               3418               0    BRCA2-201  10984
#>      start object_type         Exon                 biotype              id
#> 1 32889611  Transcript c("13", ....          protein_coding ENST00000380152
#> 2 32889642  Transcript c("GRCh3....          protein_coding ENST00000530893
#> 3 32945108  Transcript c("13", .... nonsense_mediated_decay ENST00000528762
#> 4 32953977  Transcript c("ENSE0.... nonsense_mediated_decay ENST00000470094
#> 5 32970946  Transcript c("human....         retained_intron ENST00000533776
#> 6 32889617  Transcript c("Exon"....          protein_coding ENST00000544455
#>   is_canonical                logic_name         source strand
#> 1            0 ensembl_havana_transcript ensembl_havana      1
#> 2            0    havana_homo_sapiens_37         havana      1
#> 3            0    havana_homo_sapiens_37         havana      1
#> 4            0    havana_homo_sapiens_37         havana      1
#> 5            0    havana_homo_sapiens_37         havana      1
#> 6            1   ensembl_homo_sapiens_37        ensembl      1
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $strand
#> [1] 1
#> 
#> $source
#> [1] "ensembl_havana"
#> 
```
