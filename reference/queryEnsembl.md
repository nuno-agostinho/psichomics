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
#>           gene_id canonical_transcript external_name seq_region_name version
#> 1 ENSG00000157764    ENST00000288602.6          BRAF               7       8
#> 2 ENSG00000271932    ENST00000605989.1            U6               7       1
#>   assembly_name        biotype
#> 1        GRCh37 protein_coding
#> 2        GRCh37          snRNA
#>                                                                   description
#> 1 v-raf murine sarcoma viral oncogene homolog B [Source:HGNC Symbol;Acc:1097]
#> 2                               U6 spliceosomal RNA [Source:RFAM;Acc:RF00026]
#>   feature_type     start       end              id         source
#> 1         gene 140419127 140624564 ENSG00000157764 ensembl_havana
#> 2         gene 140583872 140583978 ENSG00000271932        ensembl
#>                            logic_name strand
#> 1 ensembl_havana_gene_homo_sapiens_37     -1
#> 2               ncrna_homo_sapiens_37      1

path  <- "lookup/symbol/human/BRCA2"
query <- list(expand=1)
psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $end
#> [1] 32973805
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $db_type
#> [1] "core"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $strand
#> [1] 1
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $species
#> [1] "human"
#> 
#> $Transcript
#>        end         source strand gencode_primary db_type              id
#> 1 32973347 ensembl_havana      1               0    core ENST00000380152
#> 2 32907428         havana      1               0    core ENST00000530893
#> 3 32953632         havana      1               0    core ENST00000528762
#> 4 32972409         havana      1               0    core ENST00000470094
#> 5 32972585         havana      1               0    core ENST00000533776
#> 6 32973805        ensembl      1               0    core ENST00000544455
#>                  logic_name species seq_region_name assembly_name version
#> 1 ensembl_havana_transcript   human              13        GRCh37       3
#> 2    havana_homo_sapiens_37   human              13        GRCh37       2
#> 3    havana_homo_sapiens_37   human              13        GRCh37       1
#> 4    havana_homo_sapiens_37   human              13        GRCh37       1
#> 5    havana_homo_sapiens_37   human              13        GRCh37       1
#> 6   ensembl_homo_sapiens_37   human              13        GRCh37       1
#>      start                 biotype length         Exon is_canonical
#> 1 32889611          protein_coding  10930 c(4, 1, ....            0
#> 2 32889642          protein_coding   2011 c("human....            0
#> 3 32945108 nonsense_mediated_decay    495 c("ENSE0....            0
#> 4 32953977 nonsense_mediated_decay    842 c("13", ....            0
#> 5 32970946         retained_intron    523 c("13", ....            0
#> 6 32889617          protein_coding  10984 c("13", ....            1
#>   display_name object_type Translation.start Translation.end Translation.length
#> 1    BRCA2-001  Transcript          32890598        32972907               3418
#> 2    BRCA2-003  Transcript          32899266        32907428                481
#> 3    BRCA2-005  Transcript          32945108        32950807                 64
#> 4    BRCA2-002  Transcript          32953977        32970229                186
#> 5    BRCA2-006  Transcript                NA              NA                 NA
#> 6    BRCA2-201  Transcript          32890598        32972907               3418
#>   Translation.version Translation.species Translation.Parent
#> 1                   3               human    ENST00000380152
#> 2                   2               human    ENST00000530893
#> 3                   1               human    ENST00000528762
#> 4                   1               human    ENST00000470094
#> 5                  NA                <NA>               <NA>
#> 6                   1               human    ENST00000544455
#>   Translation.db_type Translation.object_type  Translation.id          Parent
#> 1                core             Translation ENSP00000369497 ENSG00000139618
#> 2                core             Translation ENSP00000435699 ENSG00000139618
#> 3                core             Translation ENSP00000433168 ENSG00000139618
#> 4                core             Translation ENSP00000434898 ENSG00000139618
#> 5                <NA>                    <NA>            <NA> ENSG00000139618
#> 6                core             Translation ENSP00000439902 ENSG00000139618
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $version
#> [1] 10
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $biotype
#> [1] "protein_coding"
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
```
