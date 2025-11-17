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
#>           source         gene_id feature_type
#> 1 ensembl_havana ENSG00000157764         gene
#> 2        ensembl ENSG00000271932         gene
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>       start seq_region_name        biotype external_name version strand
#> 1 140419127               7 protein_coding          BRAF       8     -1
#> 2 140583872               7          snRNA            U6       1      1
#>                id                          logic_name assembly_name
#> 1 ENSG00000157764 ensembl_havana_gene_homo_sapiens_37        GRCh37
#> 2 ENSG00000271932               ncrna_homo_sapiens_37        GRCh37
#>   canonical_transcript       end
#> 1    ENST00000288602.6 140624564
#> 2    ENST00000605989.1 140583978

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
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $Transcript
#>   object_type      end display_name species         Exon assembly_name
#> 1  Transcript 32973347    BRCA2-001   human c("GRCh3....        GRCh37
#> 2  Transcript 32907428    BRCA2-003   human c("GRCh3....        GRCh37
#> 3  Transcript 32953632    BRCA2-005   human c("GRCh3....        GRCh37
#> 4  Transcript 32972409    BRCA2-002   human c("human....        GRCh37
#> 5  Transcript 32972585    BRCA2-006   human c("GRCh3....        GRCh37
#> 6  Transcript 32973805    BRCA2-201   human c("Exon"....        GRCh37
#>   is_canonical                logic_name db_type              id strand
#> 1            0 ensembl_havana_transcript    core ENST00000380152      1
#> 2            0    havana_homo_sapiens_37    core ENST00000530893      1
#> 3            0    havana_homo_sapiens_37    core ENST00000528762      1
#> 4            0    havana_homo_sapiens_37    core ENST00000470094      1
#> 5            0    havana_homo_sapiens_37    core ENST00000533776      1
#> 6            1   ensembl_homo_sapiens_37    core ENST00000544455      1
#>   seq_region_name version    start         source length gencode_primary
#> 1              13       3 32889611 ensembl_havana  10930               0
#> 2              13       2 32889642         havana   2011               0
#> 3              13       1 32945108         havana    495               0
#> 4              13       1 32953977         havana    842               0
#> 5              13       1 32970946         havana    523               0
#> 6              13       1 32889617        ensembl  10984               0
#>            Parent Translation.species Translation.start Translation.end
#> 1 ENSG00000139618               human          32890598        32972907
#> 2 ENSG00000139618               human          32899266        32907428
#> 3 ENSG00000139618               human          32945108        32950807
#> 4 ENSG00000139618               human          32953977        32970229
#> 5 ENSG00000139618                <NA>                NA              NA
#> 6 ENSG00000139618               human          32890598        32972907
#>   Translation.object_type Translation.length  Translation.id
#> 1             Translation               3418 ENSP00000369497
#> 2             Translation                481 ENSP00000435699
#> 3             Translation                 64 ENSP00000433168
#> 4             Translation                186 ENSP00000434898
#> 5                    <NA>                 NA            <NA>
#> 6             Translation               3418 ENSP00000439902
#>   Translation.version Translation.db_type Translation.Parent
#> 1                   3                core    ENST00000380152
#> 2                   2                core    ENST00000530893
#> 3                   1                core    ENST00000528762
#> 4                   1                core    ENST00000470094
#> 5                  NA                <NA>               <NA>
#> 6                   1                core    ENST00000544455
#>                   biotype
#> 1          protein_coding
#> 2          protein_coding
#> 3 nonsense_mediated_decay
#> 4 nonsense_mediated_decay
#> 5         retained_intron
#> 6          protein_coding
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
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
#> $version
#> [1] 10
#> 
#> $seq_region_name
#> [1] "13"
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
