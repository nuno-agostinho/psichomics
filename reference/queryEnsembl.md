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
#>       start
#> 1 140419127
#> 2 140583872
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   version external_name        biotype seq_region_name
#> 1       8          BRAF protein_coding               7
#> 2       1            U6          snRNA               7
#>                            logic_name              id strand       end
#> 1 ensembl_havana_gene_homo_sapiens_37 ENSG00000157764     -1 140624564
#> 2               ncrna_homo_sapiens_37 ENSG00000271932      1 140583978
#>   canonical_transcript assembly_name         source feature_type
#> 1    ENST00000288602.6        GRCh37 ensembl_havana         gene
#> 2    ENST00000605989.1        GRCh37        ensembl         gene
#>           gene_id
#> 1 ENSG00000157764
#> 2 ENSG00000271932

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $species
#> [1] "human"
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $start
#> [1] 32889611
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
#> $end
#> [1] 32973805
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $Transcript
#>        end db_type assembly_name version length    start Translation.start
#> 1 32973347    core        GRCh37       3  10930 32889611          32890598
#> 2 32907428    core        GRCh37       2   2011 32889642          32899266
#> 3 32953632    core        GRCh37       1    495 32945108          32945108
#> 4 32972409    core        GRCh37       1    842 32953977          32953977
#> 5 32972585    core        GRCh37       1    523 32970946                NA
#> 6 32973805    core        GRCh37       1  10984 32889617          32890598
#>   Translation.length Translation.end Translation.version Translation.db_type
#> 1               3418        32972907                   3                core
#> 2                481        32907428                   2                core
#> 3                 64        32950807                   1                core
#> 4                186        32970229                   1                core
#> 5                 NA              NA                  NA                <NA>
#> 6               3418        32972907                   1                core
#>    Translation.id Translation.Parent Translation.species
#> 1 ENSP00000369497    ENST00000380152               human
#> 2 ENSP00000435699    ENST00000530893               human
#> 3 ENSP00000433168    ENST00000528762               human
#> 4 ENSP00000434898    ENST00000470094               human
#> 5            <NA>               <NA>                <NA>
#> 6 ENSP00000439902    ENST00000544455               human
#>   Translation.object_type display_name gencode_primary species seq_region_name
#> 1             Translation    BRCA2-001               0   human              13
#> 2             Translation    BRCA2-003               0   human              13
#> 3             Translation    BRCA2-005               0   human              13
#> 4             Translation    BRCA2-002               0   human              13
#> 5                    <NA>    BRCA2-006               0   human              13
#> 6             Translation    BRCA2-201               0   human              13
#>            Parent                logic_name              id is_canonical
#> 1 ENSG00000139618 ensembl_havana_transcript ENST00000380152            0
#> 2 ENSG00000139618    havana_homo_sapiens_37 ENST00000530893            0
#> 3 ENSG00000139618    havana_homo_sapiens_37 ENST00000528762            0
#> 4 ENSG00000139618    havana_homo_sapiens_37 ENST00000470094            0
#> 5 ENSG00000139618    havana_homo_sapiens_37 ENST00000533776            0
#> 6 ENSG00000139618   ensembl_homo_sapiens_37 ENST00000544455            1
#>           source strand object_type         Exon                 biotype
#> 1 ensembl_havana      1  Transcript c(328898....          protein_coding
#> 2         havana      1  Transcript c("13", ....          protein_coding
#> 3         havana      1  Transcript c("GRCh3.... nonsense_mediated_decay
#> 4         havana      1  Transcript c(1, 1, .... nonsense_mediated_decay
#> 5         havana      1  Transcript c("13", ....         retained_intron
#> 6        ensembl      1  Transcript c("13", ....          protein_coding
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $strand
#> [1] 1
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
```
