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
#>   strand                          logic_name         source              id
#> 1     -1 ensembl_havana_gene_homo_sapiens_37 ensembl_havana ENSG00000157764
#> 2      1               ncrna_homo_sapiens_37        ensembl ENSG00000271932
#>         end feature_type     start
#> 1 140624564         gene 140419127
#> 2 140583978         gene 140583872
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>          biotype assembly_name version seq_region_name canonical_transcript
#> 1 protein_coding        GRCh37       8               7    ENST00000288602.6
#> 2          snRNA        GRCh37       1               7    ENST00000605989.1
#>   external_name         gene_id
#> 1          BRAF ENSG00000157764
#> 2            U6 ENSG00000271932

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $start
#> [1] 32889611
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $end
#> [1] 32973805
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $version
#> [1] 10
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $Transcript
#>           source db_type object_type          Parent      end length
#> 1 ensembl_havana    core  Transcript ENSG00000139618 32973347  10930
#> 2         havana    core  Transcript ENSG00000139618 32907428   2011
#> 3         havana    core  Transcript ENSG00000139618 32953632    495
#> 4         havana    core  Transcript ENSG00000139618 32972409    842
#> 5         havana    core  Transcript ENSG00000139618 32972585    523
#> 6        ensembl    core  Transcript ENSG00000139618 32973805  10984
#>   display_name version    start assembly_name gencode_primary
#> 1    BRCA2-001       3 32889611        GRCh37               0
#> 2    BRCA2-003       2 32889642        GRCh37               0
#> 3    BRCA2-005       1 32945108        GRCh37               0
#> 4    BRCA2-002       1 32953977        GRCh37               0
#> 5    BRCA2-006       1 32970946        GRCh37               0
#> 6    BRCA2-201       1 32889617        GRCh37               0
#>                  logic_name                 biotype              id
#> 1 ensembl_havana_transcript          protein_coding ENST00000380152
#> 2    havana_homo_sapiens_37          protein_coding ENST00000530893
#> 3    havana_homo_sapiens_37 nonsense_mediated_decay ENST00000528762
#> 4    havana_homo_sapiens_37 nonsense_mediated_decay ENST00000470094
#> 5    havana_homo_sapiens_37         retained_intron ENST00000533776
#> 6   ensembl_homo_sapiens_37          protein_coding ENST00000544455
#>   Translation.Parent  Translation.id Translation.object_type
#> 1    ENST00000380152 ENSP00000369497             Translation
#> 2    ENST00000530893 ENSP00000435699             Translation
#> 3    ENST00000528762 ENSP00000433168             Translation
#> 4    ENST00000470094 ENSP00000434898             Translation
#> 5               <NA>            <NA>                    <NA>
#> 6    ENST00000544455 ENSP00000439902             Translation
#>   Translation.db_type Translation.species Translation.start Translation.version
#> 1                core               human          32890598                   3
#> 2                core               human          32899266                   2
#> 3                core               human          32945108                   1
#> 4                core               human          32953977                   1
#> 5                <NA>                <NA>                NA                  NA
#> 6                core               human          32890598                   1
#>   Translation.length Translation.end strand         Exon seq_region_name
#> 1               3418        32972907      1 c("13", ....              13
#> 2                481        32907428      1 c("ENSE0....              13
#> 3                 64        32950807      1 c("ENSE0....              13
#> 4                186        32970229      1 c("13", ....              13
#> 5                 NA              NA      1 c("Exon"....              13
#> 6               3418        32972907      1 c(328896....              13
#>   is_canonical species
#> 1            0   human
#> 2            0   human
#> 3            0   human
#> 4            0   human
#> 5            0   human
#> 6            1   human
#> 
#> $species
#> [1] "human"
#> 
#> $strand
#> [1] 1
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
```
