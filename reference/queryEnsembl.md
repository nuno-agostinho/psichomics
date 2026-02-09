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
#>   seq_region_name        biotype
#> 1               7 protein_coding
#> 2               7          snRNA
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   external_name         gene_id     start feature_type         source strand
#> 1          BRAF ENSG00000157764 140419127         gene ensembl_havana     -1
#> 2            U6 ENSG00000271932 140583872         gene        ensembl      1
#>                            logic_name canonical_transcript       end
#> 1 ensembl_havana_gene_homo_sapiens_37    ENST00000288602.6 140624564
#> 2               ncrna_homo_sapiens_37    ENST00000605989.1 140583978
#>                id version assembly_name
#> 1 ENSG00000157764       8        GRCh37
#> 2 ENSG00000271932       1        GRCh37

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $end
#> [1] 32973805
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $species
#> [1] "human"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $Transcript
#>                id strand seq_region_name version                logic_name
#> 1 ENST00000380152      1              13       3 ensembl_havana_transcript
#> 2 ENST00000530893      1              13       2    havana_homo_sapiens_37
#> 3 ENST00000528762      1              13       1    havana_homo_sapiens_37
#> 4 ENST00000470094      1              13       1    havana_homo_sapiens_37
#> 5 ENST00000533776      1              13       1    havana_homo_sapiens_37
#> 6 ENST00000544455      1              13       1   ensembl_homo_sapiens_37
#>   is_canonical db_type assembly_name         Exon species object_type      end
#> 1            0    core        GRCh37 c("core"....   human  Transcript 32973347
#> 2            0    core        GRCh37 c("human....   human  Transcript 32907428
#> 3            0    core        GRCh37 c("Exon"....   human  Transcript 32953632
#> 4            0    core        GRCh37 c("core"....   human  Transcript 32972409
#> 5            0    core        GRCh37 c("Exon"....   human  Transcript 32972585
#> 6            1    core        GRCh37 c("Exon"....   human  Transcript 32973805
#>   display_name Translation.version  Translation.id Translation.Parent
#> 1    BRCA2-001                   3 ENSP00000369497    ENST00000380152
#> 2    BRCA2-003                   2 ENSP00000435699    ENST00000530893
#> 3    BRCA2-005                   1 ENSP00000433168    ENST00000528762
#> 4    BRCA2-002                   1 ENSP00000434898    ENST00000470094
#> 5    BRCA2-006                  NA            <NA>               <NA>
#> 6    BRCA2-201                   1 ENSP00000439902    ENST00000544455
#>   Translation.db_type Translation.length Translation.species Translation.start
#> 1                core               3418               human          32890598
#> 2                core                481               human          32899266
#> 3                core                 64               human          32945108
#> 4                core                186               human          32953977
#> 5                <NA>                 NA                <NA>                NA
#> 6                core               3418               human          32890598
#>   Translation.end Translation.object_type                 biotype
#> 1        32972907             Translation          protein_coding
#> 2        32907428             Translation          protein_coding
#> 3        32950807             Translation nonsense_mediated_decay
#> 4        32970229             Translation nonsense_mediated_decay
#> 5              NA                    <NA>         retained_intron
#> 6        32972907             Translation          protein_coding
#>   gencode_primary          Parent length         source    start
#> 1               0 ENSG00000139618  10930 ensembl_havana 32889611
#> 2               0 ENSG00000139618   2011         havana 32889642
#> 3               0 ENSG00000139618    495         havana 32945108
#> 4               0 ENSG00000139618    842         havana 32953977
#> 5               0 ENSG00000139618    523         havana 32970946
#> 6               0 ENSG00000139618  10984        ensembl 32889617
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $db_type
#> [1] "core"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $strand
#> [1] 1
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $version
#> [1] 10
#> 
#> $start
#> [1] 32889611
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
```
