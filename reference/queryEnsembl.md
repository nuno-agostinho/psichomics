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
#>   strand              id                          logic_name         gene_id
#> 1     -1 ENSG00000157764 ensembl_havana_gene_homo_sapiens_37 ENSG00000157764
#> 2      1 ENSG00000271932               ncrna_homo_sapiens_37 ENSG00000271932
#>   seq_region_name assembly_name external_name canonical_transcript feature_type
#> 1               7        GRCh37          BRAF    ENST00000288602.6         gene
#> 2               7        GRCh37            U6    ENST00000605989.1         gene
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   version     start        biotype       end         source
#> 1       8 140419127 protein_coding 140624564 ensembl_havana
#> 2       1 140583872          snRNA 140583978        ensembl

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $display_name
#> [1] "BRCA2"
#> 
#> $strand
#> [1] 1
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $start
#> [1] 32889611
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $end
#> [1] 32973805
#> 
#> $db_type
#> [1] "core"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $version
#> [1] 10
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $Transcript
#>   length Translation.object_type Translation.Parent Translation.db_type
#> 1  10930             Translation    ENST00000380152                core
#> 2   2011             Translation    ENST00000530893                core
#> 3    495             Translation    ENST00000528762                core
#> 4    842             Translation    ENST00000470094                core
#> 5    523                    <NA>               <NA>                <NA>
#> 6  10984             Translation    ENST00000544455                core
#>   Translation.start Translation.end  Translation.id Translation.species
#> 1          32890598        32972907 ENSP00000369497               human
#> 2          32899266        32907428 ENSP00000435699               human
#> 3          32945108        32950807 ENSP00000433168               human
#> 4          32953977        32970229 ENSP00000434898               human
#> 5                NA              NA            <NA>                <NA>
#> 6          32890598        32972907 ENSP00000439902               human
#>   Translation.length Translation.version is_canonical strand display_name
#> 1               3418                   3            0      1    BRCA2-001
#> 2                481                   2            0      1    BRCA2-003
#> 3                 64                   1            0      1    BRCA2-005
#> 4                186                   1            0      1    BRCA2-002
#> 5                 NA                  NA            0      1    BRCA2-006
#> 6               3418                   1            1      1    BRCA2-201
#>   gencode_primary    start object_type seq_region_name version assembly_name
#> 1               0 32889611  Transcript              13       3        GRCh37
#> 2               0 32889642  Transcript              13       2        GRCh37
#> 3               0 32945108  Transcript              13       1        GRCh37
#> 4               0 32953977  Transcript              13       1        GRCh37
#> 5               0 32970946  Transcript              13       1        GRCh37
#> 6               0 32889617  Transcript              13       1        GRCh37
#>                   biotype         Exon species              id
#> 1          protein_coding c("13", ....   human ENST00000380152
#> 2          protein_coding c("ENSE0....   human ENST00000530893
#> 3 nonsense_mediated_decay c(1, 1, ....   human ENST00000528762
#> 4 nonsense_mediated_decay c("Exon"....   human ENST00000470094
#> 5         retained_intron c("13", ....   human ENST00000533776
#> 6          protein_coding c("ENSE0....   human ENST00000544455
#>                  logic_name         source          Parent      end db_type
#> 1 ensembl_havana_transcript ensembl_havana ENSG00000139618 32973347    core
#> 2    havana_homo_sapiens_37         havana ENSG00000139618 32907428    core
#> 3    havana_homo_sapiens_37         havana ENSG00000139618 32953632    core
#> 4    havana_homo_sapiens_37         havana ENSG00000139618 32972409    core
#> 5    havana_homo_sapiens_37         havana ENSG00000139618 32972585    core
#> 6   ensembl_homo_sapiens_37        ensembl ENSG00000139618 32973805    core
#> 
#> $species
#> [1] "human"
#> 
```
