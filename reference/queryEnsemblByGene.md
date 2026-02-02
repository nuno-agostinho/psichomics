# Query information from Ensembl

Query information from Ensembl

## Usage

``` r
queryEnsemblByGene(gene, species = NULL, assembly = NULL)

queryEnsemblByEvent(event, species = NULL, assembly = NULL, data = NULL)
```

## Arguments

- gene:

  Character: gene

- species:

  Character: species (may be `NULL` for an Ensembl identifier)

- assembly:

  Character: assembly version (may be NULL for an Ensembl identifier)

- event:

  Character: alternative splicing event

- data:

  Matrix or data frame: alternative splicing information

## Value

Information from Ensembl

## See also

Other functions to retrieve external information:
[`ensemblToUniprot()`](https://nuno-agostinho.github.io/psichomics/reference/ensemblToUniprot.md),
[`plotProtein()`](https://nuno-agostinho.github.io/psichomics/reference/plotProtein.md),
[`plotTranscripts()`](https://nuno-agostinho.github.io/psichomics/reference/plotTranscripts.md)

## Examples

``` r
queryEnsemblByGene("BRCA1", "human", "hg19")
#> $end
#> [1] 41277500
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $db_type
#> [1] "core"
#> 
#> $strand
#> [1] -1
#> 
#> $species
#> [1] "human"
#> 
#> $Transcript
#>    is_canonical species         Exon strand db_type          Parent
#> 1             0   human c(412773....     -1    core ENSG00000012048
#> 2             0   human c(412772....     -1    core ENSG00000012048
#> 3             0   human c(-1, -1....     -1    core ENSG00000012048
#> 4             0   human c("human....     -1    core ENSG00000012048
#> 5             0   human c("17", ....     -1    core ENSG00000012048
#> 6             0   human c("human....     -1    core ENSG00000012048
#> 7             1   human c("17", ....     -1    core ENSG00000012048
#> 8             0   human c(-1, -1....     -1    core ENSG00000012048
#> 9             0   human c(-1, -1....     -1    core ENSG00000012048
#> 10            0   human c(-1, -1....     -1    core ENSG00000012048
#> 11            0   human c("GRCh3....     -1    core ENSG00000012048
#> 12            0   human c(-1, -1....     -1    core ENSG00000012048
#> 13            0   human c(-1, -1....     -1    core ENSG00000012048
#> 14            0   human c(412569....     -1    core ENSG00000012048
#> 15            0   human c("17", ....     -1    core ENSG00000012048
#> 16            0   human c(-1, -1....     -1    core ENSG00000012048
#> 17            0   human c("GRCh3....     -1    core ENSG00000012048
#> 18            0   human c("core"....     -1    core ENSG00000012048
#> 19            0   human c(-1, -1....     -1    core ENSG00000012048
#> 20            0   human c(-1, -1....     -1    core ENSG00000012048
#> 21            0   human c("ENSE0....     -1    core ENSG00000012048
#> 22            0   human c("core"....     -1    core ENSG00000012048
#> 23            0   human c("core"....     -1    core ENSG00000012048
#> 24            0   human c("human....     -1    core ENSG00000012048
#> 25            0   human c(-1, -1....     -1    core ENSG00000012048
#> 26            0   human c("17", ....     -1    core ENSG00000012048
#> 27            0   human c("17", ....     -1    core ENSG00000012048
#> 28            0   human c("human....     -1    core ENSG00000012048
#> 29            0   human c(412774....     -1    core ENSG00000012048
#> 30            0   human c("Exon"....     -1    core ENSG00000012048
#> 31            0   human c("ENSE0....     -1    core ENSG00000012048
#>    seq_region_name      end    start Translation.start Translation.species
#> 1               17 41277387 41196312          41197695               human
#> 2               17 41277468 41196822          41197801               human
#> 3               17 41277346 41197580          41197695               human
#> 4               17 41277346 41197580          41197695               human
#> 5               17 41277346 41197580          41197695               human
#> 6               17 41277419 41197646          41197695               human
#> 7               17 41277500 41197646          41197695               human
#> 8               17 41277305 41197695          41256972               human
#> 9               17 41277373 41197695          41197695               human
#> 10              17 41256933 41215361          41215361               human
#> 11              17 41277376 41215361          41215361               human
#> 12              17 41277419 41215377          41215377               human
#> 13              17 41223083 41219291                NA                <NA>
#> 14              17 41256933 41228505          41228505               human
#> 15              17 41243841 41228554          41228554               human
#> 16              17 41277332 41243115                NA                <NA>
#> 17              17 41247883 41245587          41245587               human
#> 18              17 41277340 41245601          41245601               human
#> 19              17 41277381 41245603          41245603               human
#> 20              17 41277317 41246129          41262552               human
#> 21              17 41277317 41246129          41246129               human
#> 22              17 41277467 41246129          41246129               human
#> 23              17 41256908 41246187          41246187               human
#> 24              17 41277370 41247863          41247863               human
#> 25              17 41277387 41251848          41256972               human
#> 26              17 41277338 41256206          41256206               human
#> 27              17 41277500 41196313          41197695               human
#> 28              17 41277500 41196313          41197695               human
#> 29              17 41277468 41196313          41197695               human
#> 30              17 41276132 41196313          41197695               human
#> 31              17 41277468 41196313          41197695               human
#>    Translation.length Translation.version Translation.db_type
#> 1                1863                   3                core
#> 2                 699                   1                core
#> 3                 173                   1                core
#> 4                 354                   1                core
#> 5                  96                   1                core
#> 6                1816                   1                core
#> 7                1884                   2                core
#> 8                  63                   1                core
#> 9                 759                   2                core
#> 10                498                   1                core
#> 11                623                   1                core
#> 12                572                   1                core
#> 13                 NA                  NA                <NA>
#> 14                266                   1                core
#> 15                242                   1                core
#> 16                 NA                  NA                <NA>
#> 17                437                   3                core
#> 18                649                   1                core
#> 19                622                   1                core
#> 20                 59                   1                core
#> 21                177                   1                core
#> 22                473                   1                core
#> 23                319                   1                core
#> 24                222                   1                core
#> 25                 63                   1                core
#> 26                 98                   1                core
#> 27               1598                   6                core
#> 28                721                   5                core
#> 29               1624                   4                core
#> 30                680                   3                core
#> 31               1567                   4                core
#>    Translation.object_type Translation.Parent  Translation.id Translation.end
#> 1              Translation    ENST00000357654 ENSP00000350283        41276113
#> 2              Translation    ENST00000468300 ENSP00000417148        41276113
#> 3              Translation    ENST00000586385 ENSP00000465818        41277202
#> 4              Translation    ENST00000591534 ENSP00000467329        41226495
#> 5              Translation    ENST00000591849 ENSP00000465347        41202109
#> 6              Translation    ENST00000493795 ENSP00000418775        41258543
#> 7              Translation    ENST00000471181 ENSP00000418960        41276113
#> 8              Translation    ENST00000461221 ENSP00000418548        41276113
#> 9              Translation    ENST00000491747 ENSP00000420705        41276113
#> 10             Translation    ENST00000484087 ENSP00000419481        41256933
#> 11             Translation    ENST00000478531 ENSP00000420412        41276113
#> 12             Translation    ENST00000493919 ENSP00000418819        41258543
#> 13                    <NA>               <NA>            <NA>              NA
#> 14             Translation    ENST00000487825 ENSP00000418212        41256933
#> 15             Translation    ENST00000461574 ENSP00000417241        41243841
#> 16                    <NA>               <NA>            <NA>              NA
#> 17             Translation    ENST00000412061 ENSP00000397145        41247883
#> 18             Translation    ENST00000470026 ENSP00000419274        41276113
#> 19             Translation    ENST00000477152 ENSP00000419988        41276113
#> 20             Translation    ENST00000492859 ENSP00000420253        41276113
#> 21             Translation    ENST00000497488 ENSP00000418986        41246659
#> 22             Translation    ENST00000494123 ENSP00000419103        41276113
#> 23             Translation    ENST00000473961 ENSP00000420201        41256908
#> 24             Translation    ENST00000476777 ENSP00000417554        41276113
#> 25             Translation    ENST00000461798 ENSP00000417988        41276113
#> 26             Translation    ENST00000489037 ENSP00000420781        41276113
#> 27             Translation    ENST00000354071 ENSP00000326002        41276113
#> 28             Translation    ENST00000352993 ENSP00000312236        41276113
#> 29             Translation    ENST00000346315 ENSP00000246907        41276113
#> 30             Translation    ENST00000351666 ENSP00000338007        41276113
#> 31             Translation    ENST00000309486 ENSP00000310938        41246659
#>    display_name length         source assembly_name version
#> 1     BRCA1-001   7094 ensembl_havana        GRCh37       3
#> 2     BRCA1-007   3273 ensembl_havana        GRCh37       1
#> 3     BRCA1-023    781         havana        GRCh37       1
#> 4     BRCA1-024   1282         havana        GRCh37       1
#> 5     BRCA1-025    563         havana        GRCh37       1
#> 6     BRCA1-006   5732 ensembl_havana        GRCh37       1
#> 7     BRCA1-005   5936 ensembl_havana        GRCh37       2
#> 8     BRCA1-010   5693         havana        GRCh37       1
#> 9     BRCA1-014   2379         havana        GRCh37       2
#> 10    BRCA1-015   1495         havana        GRCh37       1
#> 11    BRCA1-009   1972         havana        GRCh37       1
#> 12    BRCA1-008   1948         havana        GRCh37       1
#> 13    BRCA1-021    561         havana        GRCh37       1
#> 14    BRCA1-019    800         havana        GRCh37       1
#> 15    BRCA1-022    726         havana        GRCh37       1
#> 16    BRCA1-012   4497         havana        GRCh37       1
#> 17    BRCA1-026   1312         havana        GRCh37       3
#> 18    BRCA1-011   2108         havana        GRCh37       1
#> 19    BRCA1-004   1980         havana        GRCh37       1
#> 20    BRCA1-002   1584         havana        GRCh37       1
#> 21    BRCA1-003    779         havana        GRCh37       1
#> 22    BRCA1-013   1612         havana        GRCh37       1
#> 23    BRCA1-018    958         havana        GRCh37       1
#> 24    BRCA1-017    769         havana        GRCh37       1
#> 25    BRCA1-020    582         havana        GRCh37       1
#> 26    BRCA1-016    455         havana        GRCh37       1
#> 27    BRCA1-205   6411        ensembl        GRCh37       3
#> 28    BRCA1-204   3780        ensembl        GRCh37       3
#> 29    BRCA1-202   6451        ensembl        GRCh37       3
#> 30    BRCA1-203   3444        ensembl        GRCh37       3
#> 31    BRCA1-201   7114        ensembl        GRCh37       4
#>                    biotype object_type                logic_name
#> 1           protein_coding  Transcript ensembl_havana_transcript
#> 2           protein_coding  Transcript ensembl_havana_transcript
#> 3           protein_coding  Transcript    havana_homo_sapiens_37
#> 4           protein_coding  Transcript    havana_homo_sapiens_37
#> 5           protein_coding  Transcript    havana_homo_sapiens_37
#> 6           protein_coding  Transcript ensembl_havana_transcript
#> 7           protein_coding  Transcript ensembl_havana_transcript
#> 8  nonsense_mediated_decay  Transcript    havana_homo_sapiens_37
#> 9           protein_coding  Transcript    havana_homo_sapiens_37
#> 10          protein_coding  Transcript    havana_homo_sapiens_37
#> 11          protein_coding  Transcript    havana_homo_sapiens_37
#> 12          protein_coding  Transcript    havana_homo_sapiens_37
#> 13         retained_intron  Transcript    havana_homo_sapiens_37
#> 14          protein_coding  Transcript    havana_homo_sapiens_37
#> 15          protein_coding  Transcript    havana_homo_sapiens_37
#> 16         retained_intron  Transcript    havana_homo_sapiens_37
#> 17          non_stop_decay  Transcript    havana_homo_sapiens_37
#> 18          protein_coding  Transcript    havana_homo_sapiens_37
#> 19          protein_coding  Transcript    havana_homo_sapiens_37
#> 20 nonsense_mediated_decay  Transcript    havana_homo_sapiens_37
#> 21          protein_coding  Transcript    havana_homo_sapiens_37
#> 22          protein_coding  Transcript    havana_homo_sapiens_37
#> 23          protein_coding  Transcript    havana_homo_sapiens_37
#> 24          protein_coding  Transcript    havana_homo_sapiens_37
#> 25 nonsense_mediated_decay  Transcript    havana_homo_sapiens_37
#> 26          protein_coding  Transcript    havana_homo_sapiens_37
#> 27          protein_coding  Transcript   ensembl_homo_sapiens_37
#> 28          protein_coding  Transcript   ensembl_homo_sapiens_37
#> 29          protein_coding  Transcript   ensembl_homo_sapiens_37
#> 30          protein_coding  Transcript   ensembl_homo_sapiens_37
#> 31          protein_coding  Transcript   ensembl_homo_sapiens_37
#>    gencode_primary              id
#> 1                0 ENST00000357654
#> 2                0 ENST00000468300
#> 3                0 ENST00000586385
#> 4                0 ENST00000591534
#> 5                0 ENST00000591849
#> 6                0 ENST00000493795
#> 7                0 ENST00000471181
#> 8                0 ENST00000461221
#> 9                0 ENST00000491747
#> 10               0 ENST00000484087
#> 11               0 ENST00000478531
#> 12               0 ENST00000493919
#> 13               0 ENST00000472490
#> 14               0 ENST00000487825
#> 15               0 ENST00000461574
#> 16               0 ENST00000467274
#> 17               0 ENST00000412061
#> 18               0 ENST00000470026
#> 19               0 ENST00000477152
#> 20               0 ENST00000492859
#> 21               0 ENST00000497488
#> 22               0 ENST00000494123
#> 23               0 ENST00000473961
#> 24               0 ENST00000476777
#> 25               0 ENST00000461798
#> 26               0 ENST00000489037
#> 27               0 ENST00000354071
#> 28               0 ENST00000352993
#> 29               0 ENST00000346315
#> 30               0 ENST00000351666
#> 31               0 ENST00000309486
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
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
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $version
#> [1] 15
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $start
#> [1] 41196312
#> 
queryEnsemblByGene("ENSG00000139618")
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $strand
#> [1] 1
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $end
#> [1] 32973805
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $start
#> [1] 32889611
#> 
#> $Transcript
#>   display_name assembly_name              id version gencode_primary    start
#> 1    BRCA2-001        GRCh37 ENST00000380152       3               0 32889611
#> 2    BRCA2-003        GRCh37 ENST00000530893       2               0 32889642
#> 3    BRCA2-005        GRCh37 ENST00000528762       1               0 32945108
#> 4    BRCA2-002        GRCh37 ENST00000470094       1               0 32953977
#> 5    BRCA2-006        GRCh37 ENST00000533776       1               0 32970946
#> 6    BRCA2-201        GRCh37 ENST00000544455       1               0 32889617
#>   length db_type  Translation.id Translation.end Translation.Parent
#> 1  10930    core ENSP00000369497        32972907    ENST00000380152
#> 2   2011    core ENSP00000435699        32907428    ENST00000530893
#> 3    495    core ENSP00000433168        32950807    ENST00000528762
#> 4    842    core ENSP00000434898        32970229    ENST00000470094
#> 5    523    core            <NA>              NA               <NA>
#> 6  10984    core ENSP00000439902        32972907    ENST00000544455
#>   Translation.start Translation.length Translation.version Translation.species
#> 1          32890598               3418                   3        homo_sapiens
#> 2          32899266                481                   2        homo_sapiens
#> 3          32945108                 64                   1        homo_sapiens
#> 4          32953977                186                   1        homo_sapiens
#> 5                NA                 NA                  NA                <NA>
#> 6          32890598               3418                   1        homo_sapiens
#>   Translation.db_type Translation.object_type object_type         Exon      end
#> 1                core             Translation  Transcript c(1, 1, .... 32973347
#> 2                core             Translation  Transcript c("Exon".... 32907428
#> 3                core             Translation  Transcript c(329451.... 32953632
#> 4                core             Translation  Transcript c(329540.... 32972409
#> 5                <NA>                    <NA>  Transcript c("GRCh3.... 32972585
#> 6                core             Translation  Transcript c("ENSE0.... 32973805
#>                  logic_name strand      species                 biotype
#> 1 ensembl_havana_transcript      1 homo_sapiens          protein_coding
#> 2    havana_homo_sapiens_37      1 homo_sapiens          protein_coding
#> 3    havana_homo_sapiens_37      1 homo_sapiens nonsense_mediated_decay
#> 4    havana_homo_sapiens_37      1 homo_sapiens nonsense_mediated_decay
#> 5    havana_homo_sapiens_37      1 homo_sapiens         retained_intron
#> 6   ensembl_homo_sapiens_37      1 homo_sapiens          protein_coding
#>   seq_region_name          Parent is_canonical         source
#> 1              13 ENSG00000139618            0 ensembl_havana
#> 2              13 ENSG00000139618            0         havana
#> 3              13 ENSG00000139618            0         havana
#> 4              13 ENSG00000139618            0         havana
#> 5              13 ENSG00000139618            0         havana
#> 6              13 ENSG00000139618            1        ensembl
#> 
#> $version
#> [1] 10
#> 
#> $db_type
#> [1] "core"
#> 
#> $object_type
#> [1] "Gene"
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $db_type
#> [1] "core"
#> 
#> $version
#> [1] 15
#> 
#> $end
#> [1] 41277500
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $start
#> [1] 41196312
#> 
#> $species
#> [1] "human"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $strand
#> [1] -1
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $Transcript
#>    species          Parent seq_region_name assembly_name version db_type
#> 1    human ENSG00000012048              17        GRCh37       3    core
#> 2    human ENSG00000012048              17        GRCh37       1    core
#> 3    human ENSG00000012048              17        GRCh37       1    core
#> 4    human ENSG00000012048              17        GRCh37       1    core
#> 5    human ENSG00000012048              17        GRCh37       1    core
#> 6    human ENSG00000012048              17        GRCh37       1    core
#> 7    human ENSG00000012048              17        GRCh37       2    core
#> 8    human ENSG00000012048              17        GRCh37       1    core
#> 9    human ENSG00000012048              17        GRCh37       2    core
#> 10   human ENSG00000012048              17        GRCh37       1    core
#> 11   human ENSG00000012048              17        GRCh37       1    core
#> 12   human ENSG00000012048              17        GRCh37       1    core
#> 13   human ENSG00000012048              17        GRCh37       1    core
#> 14   human ENSG00000012048              17        GRCh37       1    core
#> 15   human ENSG00000012048              17        GRCh37       1    core
#> 16   human ENSG00000012048              17        GRCh37       1    core
#> 17   human ENSG00000012048              17        GRCh37       3    core
#> 18   human ENSG00000012048              17        GRCh37       1    core
#> 19   human ENSG00000012048              17        GRCh37       1    core
#> 20   human ENSG00000012048              17        GRCh37       1    core
#> 21   human ENSG00000012048              17        GRCh37       1    core
#> 22   human ENSG00000012048              17        GRCh37       1    core
#> 23   human ENSG00000012048              17        GRCh37       1    core
#> 24   human ENSG00000012048              17        GRCh37       1    core
#> 25   human ENSG00000012048              17        GRCh37       1    core
#> 26   human ENSG00000012048              17        GRCh37       1    core
#> 27   human ENSG00000012048              17        GRCh37       3    core
#> 28   human ENSG00000012048              17        GRCh37       3    core
#> 29   human ENSG00000012048              17        GRCh37       3    core
#> 30   human ENSG00000012048              17        GRCh37       3    core
#> 31   human ENSG00000012048              17        GRCh37       4    core
#>         end display_name gencode_primary Translation.start Translation.length
#> 1  41277387    BRCA1-001               0          41197695               1863
#> 2  41277468    BRCA1-007               0          41197801                699
#> 3  41277346    BRCA1-023               0          41197695                173
#> 4  41277346    BRCA1-024               0          41197695                354
#> 5  41277346    BRCA1-025               0          41197695                 96
#> 6  41277419    BRCA1-006               0          41197695               1816
#> 7  41277500    BRCA1-005               0          41197695               1884
#> 8  41277305    BRCA1-010               0          41256972                 63
#> 9  41277373    BRCA1-014               0          41197695                759
#> 10 41256933    BRCA1-015               0          41215361                498
#> 11 41277376    BRCA1-009               0          41215361                623
#> 12 41277419    BRCA1-008               0          41215377                572
#> 13 41223083    BRCA1-021               0                NA                 NA
#> 14 41256933    BRCA1-019               0          41228505                266
#> 15 41243841    BRCA1-022               0          41228554                242
#> 16 41277332    BRCA1-012               0                NA                 NA
#> 17 41247883    BRCA1-026               0          41245587                437
#> 18 41277340    BRCA1-011               0          41245601                649
#> 19 41277381    BRCA1-004               0          41245603                622
#> 20 41277317    BRCA1-002               0          41262552                 59
#> 21 41277317    BRCA1-003               0          41246129                177
#> 22 41277467    BRCA1-013               0          41246129                473
#> 23 41256908    BRCA1-018               0          41246187                319
#> 24 41277370    BRCA1-017               0          41247863                222
#> 25 41277387    BRCA1-020               0          41256972                 63
#> 26 41277338    BRCA1-016               0          41256206                 98
#> 27 41277500    BRCA1-205               0          41197695               1598
#> 28 41277500    BRCA1-204               0          41197695                721
#> 29 41277468    BRCA1-202               0          41197695               1624
#> 30 41276132    BRCA1-203               0          41197695                680
#> 31 41277468    BRCA1-201               0          41197695               1567
#>     Translation.id Translation.version Translation.db_type Translation.end
#> 1  ENSP00000350283                   3                core        41276113
#> 2  ENSP00000417148                   1                core        41276113
#> 3  ENSP00000465818                   1                core        41277202
#> 4  ENSP00000467329                   1                core        41226495
#> 5  ENSP00000465347                   1                core        41202109
#> 6  ENSP00000418775                   1                core        41258543
#> 7  ENSP00000418960                   2                core        41276113
#> 8  ENSP00000418548                   1                core        41276113
#> 9  ENSP00000420705                   2                core        41276113
#> 10 ENSP00000419481                   1                core        41256933
#> 11 ENSP00000420412                   1                core        41276113
#> 12 ENSP00000418819                   1                core        41258543
#> 13            <NA>                  NA                <NA>              NA
#> 14 ENSP00000418212                   1                core        41256933
#> 15 ENSP00000417241                   1                core        41243841
#> 16            <NA>                  NA                <NA>              NA
#> 17 ENSP00000397145                   3                core        41247883
#> 18 ENSP00000419274                   1                core        41276113
#> 19 ENSP00000419988                   1                core        41276113
#> 20 ENSP00000420253                   1                core        41276113
#> 21 ENSP00000418986                   1                core        41246659
#> 22 ENSP00000419103                   1                core        41276113
#> 23 ENSP00000420201                   1                core        41256908
#> 24 ENSP00000417554                   1                core        41276113
#> 25 ENSP00000417988                   1                core        41276113
#> 26 ENSP00000420781                   1                core        41276113
#> 27 ENSP00000326002                   6                core        41276113
#> 28 ENSP00000312236                   5                core        41276113
#> 29 ENSP00000246907                   4                core        41276113
#> 30 ENSP00000338007                   3                core        41276113
#> 31 ENSP00000310938                   4                core        41246659
#>    Translation.Parent Translation.object_type Translation.species length
#> 1     ENST00000357654             Translation               human   7094
#> 2     ENST00000468300             Translation               human   3273
#> 3     ENST00000586385             Translation               human    781
#> 4     ENST00000591534             Translation               human   1282
#> 5     ENST00000591849             Translation               human    563
#> 6     ENST00000493795             Translation               human   5732
#> 7     ENST00000471181             Translation               human   5936
#> 8     ENST00000461221             Translation               human   5693
#> 9     ENST00000491747             Translation               human   2379
#> 10    ENST00000484087             Translation               human   1495
#> 11    ENST00000478531             Translation               human   1972
#> 12    ENST00000493919             Translation               human   1948
#> 13               <NA>                    <NA>                <NA>    561
#> 14    ENST00000487825             Translation               human    800
#> 15    ENST00000461574             Translation               human    726
#> 16               <NA>                    <NA>                <NA>   4497
#> 17    ENST00000412061             Translation               human   1312
#> 18    ENST00000470026             Translation               human   2108
#> 19    ENST00000477152             Translation               human   1980
#> 20    ENST00000492859             Translation               human   1584
#> 21    ENST00000497488             Translation               human    779
#> 22    ENST00000494123             Translation               human   1612
#> 23    ENST00000473961             Translation               human    958
#> 24    ENST00000476777             Translation               human    769
#> 25    ENST00000461798             Translation               human    582
#> 26    ENST00000489037             Translation               human    455
#> 27    ENST00000354071             Translation               human   6411
#> 28    ENST00000352993             Translation               human   3780
#> 29    ENST00000346315             Translation               human   6451
#> 30    ENST00000351666             Translation               human   3444
#> 31    ENST00000309486             Translation               human   7114
#>       start         Exon object_type                 biotype is_canonical
#> 1  41196312 c("Exon"....  Transcript          protein_coding            0
#> 2  41196822 c("17", ....  Transcript          protein_coding            0
#> 3  41197580 c("Exon"....  Transcript          protein_coding            0
#> 4  41197580 c(-1, -1....  Transcript          protein_coding            0
#> 5  41197580 c("human....  Transcript          protein_coding            0
#> 6  41197646 c("human....  Transcript          protein_coding            0
#> 7  41197646 c("Exon"....  Transcript          protein_coding            1
#> 8  41197695 c(412773....  Transcript nonsense_mediated_decay            0
#> 9  41197695 c("17", ....  Transcript          protein_coding            0
#> 10 41215361 c("17", ....  Transcript          protein_coding            0
#> 11 41215361 c("17", ....  Transcript          protein_coding            0
#> 12 41215377 c(412774....  Transcript          protein_coding            0
#> 13 41219291 c("ENSE0....  Transcript         retained_intron            0
#> 14 41228505 c(412569....  Transcript          protein_coding            0
#> 15 41228554 c(-1, -1....  Transcript          protein_coding            0
#> 16 41243115 c(412772....  Transcript         retained_intron            0
#> 17 41245587 c("human....  Transcript          non_stop_decay            0
#> 18 41245601 c("17", ....  Transcript          protein_coding            0
#> 19 41245603 c(412772....  Transcript          protein_coding            0
#> 20 41246129 c(412772....  Transcript nonsense_mediated_decay            0
#> 21 41246129 c("ENSE0....  Transcript          protein_coding            0
#> 22 41246129 c("17", ....  Transcript          protein_coding            0
#> 23 41246187 c(412569....  Transcript          protein_coding            0
#> 24 41247863 c("17", ....  Transcript          protein_coding            0
#> 25 41251848 c(412773....  Transcript nonsense_mediated_decay            0
#> 26 41256206 c("ENSE0....  Transcript          protein_coding            0
#> 27 41196313 c(412775....  Transcript          protein_coding            0
#> 28 41196313 c("human....  Transcript          protein_coding            0
#> 29 41196313 c(-1, -1....  Transcript          protein_coding            0
#> 30 41196313 c("Exon"....  Transcript          protein_coding            0
#> 31 41196313 c("17", ....  Transcript          protein_coding            0
#>                 id                logic_name strand         source
#> 1  ENST00000357654 ensembl_havana_transcript     -1 ensembl_havana
#> 2  ENST00000468300 ensembl_havana_transcript     -1 ensembl_havana
#> 3  ENST00000586385    havana_homo_sapiens_37     -1         havana
#> 4  ENST00000591534    havana_homo_sapiens_37     -1         havana
#> 5  ENST00000591849    havana_homo_sapiens_37     -1         havana
#> 6  ENST00000493795 ensembl_havana_transcript     -1 ensembl_havana
#> 7  ENST00000471181 ensembl_havana_transcript     -1 ensembl_havana
#> 8  ENST00000461221    havana_homo_sapiens_37     -1         havana
#> 9  ENST00000491747    havana_homo_sapiens_37     -1         havana
#> 10 ENST00000484087    havana_homo_sapiens_37     -1         havana
#> 11 ENST00000478531    havana_homo_sapiens_37     -1         havana
#> 12 ENST00000493919    havana_homo_sapiens_37     -1         havana
#> 13 ENST00000472490    havana_homo_sapiens_37     -1         havana
#> 14 ENST00000487825    havana_homo_sapiens_37     -1         havana
#> 15 ENST00000461574    havana_homo_sapiens_37     -1         havana
#> 16 ENST00000467274    havana_homo_sapiens_37     -1         havana
#> 17 ENST00000412061    havana_homo_sapiens_37     -1         havana
#> 18 ENST00000470026    havana_homo_sapiens_37     -1         havana
#> 19 ENST00000477152    havana_homo_sapiens_37     -1         havana
#> 20 ENST00000492859    havana_homo_sapiens_37     -1         havana
#> 21 ENST00000497488    havana_homo_sapiens_37     -1         havana
#> 22 ENST00000494123    havana_homo_sapiens_37     -1         havana
#> 23 ENST00000473961    havana_homo_sapiens_37     -1         havana
#> 24 ENST00000476777    havana_homo_sapiens_37     -1         havana
#> 25 ENST00000461798    havana_homo_sapiens_37     -1         havana
#> 26 ENST00000489037    havana_homo_sapiens_37     -1         havana
#> 27 ENST00000354071   ensembl_homo_sapiens_37     -1        ensembl
#> 28 ENST00000352993   ensembl_homo_sapiens_37     -1        ensembl
#> 29 ENST00000346315   ensembl_homo_sapiens_37     -1        ensembl
#> 30 ENST00000351666   ensembl_homo_sapiens_37     -1        ensembl
#> 31 ENST00000309486   ensembl_homo_sapiens_37     -1        ensembl
#> 
```
