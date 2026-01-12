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
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $version
#> [1] 15
#> 
#> $start
#> [1] 41196312
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $end
#> [1] 41277500
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $strand
#> [1] -1
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $db_type
#> [1] "core"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $species
#> [1] "human"
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $Transcript
#>    Translation.version Translation.start Translation.length Translation.end
#> 1                    3          41197695               1863        41276113
#> 2                    1          41197801                699        41276113
#> 3                    1          41197695                173        41277202
#> 4                    1          41197695                354        41226495
#> 5                    1          41197695                 96        41202109
#> 6                    1          41197695               1816        41258543
#> 7                    2          41197695               1884        41276113
#> 8                    1          41256972                 63        41276113
#> 9                    2          41197695                759        41276113
#> 10                   1          41215361                498        41256933
#> 11                   1          41215361                623        41276113
#> 12                   1          41215377                572        41258543
#> 13                  NA                NA                 NA              NA
#> 14                   1          41228505                266        41256933
#> 15                   1          41228554                242        41243841
#> 16                  NA                NA                 NA              NA
#> 17                   3          41245587                437        41247883
#> 18                   1          41245601                649        41276113
#> 19                   1          41245603                622        41276113
#> 20                   1          41262552                 59        41276113
#> 21                   1          41246129                177        41246659
#> 22                   1          41246129                473        41276113
#> 23                   1          41246187                319        41256908
#> 24                   1          41247863                222        41276113
#> 25                   1          41256972                 63        41276113
#> 26                   1          41256206                 98        41276113
#> 27                   6          41197695               1598        41276113
#> 28                   5          41197695                721        41276113
#> 29                   4          41197695               1624        41276113
#> 30                   3          41197695                680        41276113
#> 31                   4          41197695               1567        41246659
#>    Translation.object_type Translation.db_type  Translation.id
#> 1              Translation                core ENSP00000350283
#> 2              Translation                core ENSP00000417148
#> 3              Translation                core ENSP00000465818
#> 4              Translation                core ENSP00000467329
#> 5              Translation                core ENSP00000465347
#> 6              Translation                core ENSP00000418775
#> 7              Translation                core ENSP00000418960
#> 8              Translation                core ENSP00000418548
#> 9              Translation                core ENSP00000420705
#> 10             Translation                core ENSP00000419481
#> 11             Translation                core ENSP00000420412
#> 12             Translation                core ENSP00000418819
#> 13                    <NA>                <NA>            <NA>
#> 14             Translation                core ENSP00000418212
#> 15             Translation                core ENSP00000417241
#> 16                    <NA>                <NA>            <NA>
#> 17             Translation                core ENSP00000397145
#> 18             Translation                core ENSP00000419274
#> 19             Translation                core ENSP00000419988
#> 20             Translation                core ENSP00000420253
#> 21             Translation                core ENSP00000418986
#> 22             Translation                core ENSP00000419103
#> 23             Translation                core ENSP00000420201
#> 24             Translation                core ENSP00000417554
#> 25             Translation                core ENSP00000417988
#> 26             Translation                core ENSP00000420781
#> 27             Translation                core ENSP00000326002
#> 28             Translation                core ENSP00000312236
#> 29             Translation                core ENSP00000246907
#> 30             Translation                core ENSP00000338007
#> 31             Translation                core ENSP00000310938
#>    Translation.species Translation.Parent display_name object_type is_canonical
#> 1                human    ENST00000357654    BRCA1-001  Transcript            0
#> 2                human    ENST00000468300    BRCA1-007  Transcript            0
#> 3                human    ENST00000586385    BRCA1-023  Transcript            0
#> 4                human    ENST00000591534    BRCA1-024  Transcript            0
#> 5                human    ENST00000591849    BRCA1-025  Transcript            0
#> 6                human    ENST00000493795    BRCA1-006  Transcript            0
#> 7                human    ENST00000471181    BRCA1-005  Transcript            1
#> 8                human    ENST00000461221    BRCA1-010  Transcript            0
#> 9                human    ENST00000491747    BRCA1-014  Transcript            0
#> 10               human    ENST00000484087    BRCA1-015  Transcript            0
#> 11               human    ENST00000478531    BRCA1-009  Transcript            0
#> 12               human    ENST00000493919    BRCA1-008  Transcript            0
#> 13                <NA>               <NA>    BRCA1-021  Transcript            0
#> 14               human    ENST00000487825    BRCA1-019  Transcript            0
#> 15               human    ENST00000461574    BRCA1-022  Transcript            0
#> 16                <NA>               <NA>    BRCA1-012  Transcript            0
#> 17               human    ENST00000412061    BRCA1-026  Transcript            0
#> 18               human    ENST00000470026    BRCA1-011  Transcript            0
#> 19               human    ENST00000477152    BRCA1-004  Transcript            0
#> 20               human    ENST00000492859    BRCA1-002  Transcript            0
#> 21               human    ENST00000497488    BRCA1-003  Transcript            0
#> 22               human    ENST00000494123    BRCA1-013  Transcript            0
#> 23               human    ENST00000473961    BRCA1-018  Transcript            0
#> 24               human    ENST00000476777    BRCA1-017  Transcript            0
#> 25               human    ENST00000461798    BRCA1-020  Transcript            0
#> 26               human    ENST00000489037    BRCA1-016  Transcript            0
#> 27               human    ENST00000354071    BRCA1-205  Transcript            0
#> 28               human    ENST00000352993    BRCA1-204  Transcript            0
#> 29               human    ENST00000346315    BRCA1-202  Transcript            0
#> 30               human    ENST00000351666    BRCA1-203  Transcript            0
#> 31               human    ENST00000309486    BRCA1-201  Transcript            0
#>            Exon          Parent version length                 biotype    start
#> 1  c(1, 1, .... ENSG00000012048       3   7094          protein_coding 41196312
#> 2  c("human.... ENSG00000012048       1   3273          protein_coding 41196822
#> 3  c("Exon".... ENSG00000012048       1    781          protein_coding 41197580
#> 4  c("17", .... ENSG00000012048       1   1282          protein_coding 41197580
#> 5  c("GRCh3.... ENSG00000012048       1    563          protein_coding 41197580
#> 6  c(1, 1, .... ENSG00000012048       1   5732          protein_coding 41197646
#> 7  c("17", .... ENSG00000012048       2   5936          protein_coding 41197646
#> 8  c(412773.... ENSG00000012048       1   5693 nonsense_mediated_decay 41197695
#> 9  c(412773.... ENSG00000012048       2   2379          protein_coding 41197695
#> 10 c(1, 1, .... ENSG00000012048       1   1495          protein_coding 41215361
#> 11 c(1, 1, .... ENSG00000012048       1   1972          protein_coding 41215361
#> 12 c("Exon".... ENSG00000012048       1   1948          protein_coding 41215377
#> 13 c("GRCh3.... ENSG00000012048       1    561         retained_intron 41219291
#> 14 c(1, 1, .... ENSG00000012048       1    800          protein_coding 41228505
#> 15 c(412438.... ENSG00000012048       1    726          protein_coding 41228554
#> 16 c(412772.... ENSG00000012048       1   4497         retained_intron 41243115
#> 17 c(412478.... ENSG00000012048       3   1312          non_stop_decay 41245587
#> 18 c("core".... ENSG00000012048       1   2108          protein_coding 41245601
#> 19 c("human.... ENSG00000012048       1   1980          protein_coding 41245603
#> 20 c("Exon".... ENSG00000012048       1   1584 nonsense_mediated_decay 41246129
#> 21 c(1, 1),.... ENSG00000012048       1    779          protein_coding 41246129
#> 22 c(1, 1, .... ENSG00000012048       1   1612          protein_coding 41246129
#> 23 c("17", .... ENSG00000012048       1    958          protein_coding 41246187
#> 24 c(-1, -1.... ENSG00000012048       1    769          protein_coding 41247863
#> 25 c(1, 1, .... ENSG00000012048       1    582 nonsense_mediated_decay 41251848
#> 26 c(1, 1, .... ENSG00000012048       1    455          protein_coding 41256206
#> 27 c(-1, -1.... ENSG00000012048       3   6411          protein_coding 41196313
#> 28 c("17", .... ENSG00000012048       3   3780          protein_coding 41196313
#> 29 c(1, 1, .... ENSG00000012048       3   6451          protein_coding 41196313
#> 30 c("17", .... ENSG00000012048       3   3444          protein_coding 41196313
#> 31 c(-1, -1.... ENSG00000012048       4   7114          protein_coding 41196313
#>                   logic_name              id gencode_primary db_type strand
#> 1  ensembl_havana_transcript ENST00000357654               0    core     -1
#> 2  ensembl_havana_transcript ENST00000468300               0    core     -1
#> 3     havana_homo_sapiens_37 ENST00000586385               0    core     -1
#> 4     havana_homo_sapiens_37 ENST00000591534               0    core     -1
#> 5     havana_homo_sapiens_37 ENST00000591849               0    core     -1
#> 6  ensembl_havana_transcript ENST00000493795               0    core     -1
#> 7  ensembl_havana_transcript ENST00000471181               0    core     -1
#> 8     havana_homo_sapiens_37 ENST00000461221               0    core     -1
#> 9     havana_homo_sapiens_37 ENST00000491747               0    core     -1
#> 10    havana_homo_sapiens_37 ENST00000484087               0    core     -1
#> 11    havana_homo_sapiens_37 ENST00000478531               0    core     -1
#> 12    havana_homo_sapiens_37 ENST00000493919               0    core     -1
#> 13    havana_homo_sapiens_37 ENST00000472490               0    core     -1
#> 14    havana_homo_sapiens_37 ENST00000487825               0    core     -1
#> 15    havana_homo_sapiens_37 ENST00000461574               0    core     -1
#> 16    havana_homo_sapiens_37 ENST00000467274               0    core     -1
#> 17    havana_homo_sapiens_37 ENST00000412061               0    core     -1
#> 18    havana_homo_sapiens_37 ENST00000470026               0    core     -1
#> 19    havana_homo_sapiens_37 ENST00000477152               0    core     -1
#> 20    havana_homo_sapiens_37 ENST00000492859               0    core     -1
#> 21    havana_homo_sapiens_37 ENST00000497488               0    core     -1
#> 22    havana_homo_sapiens_37 ENST00000494123               0    core     -1
#> 23    havana_homo_sapiens_37 ENST00000473961               0    core     -1
#> 24    havana_homo_sapiens_37 ENST00000476777               0    core     -1
#> 25    havana_homo_sapiens_37 ENST00000461798               0    core     -1
#> 26    havana_homo_sapiens_37 ENST00000489037               0    core     -1
#> 27   ensembl_homo_sapiens_37 ENST00000354071               0    core     -1
#> 28   ensembl_homo_sapiens_37 ENST00000352993               0    core     -1
#> 29   ensembl_homo_sapiens_37 ENST00000346315               0    core     -1
#> 30   ensembl_homo_sapiens_37 ENST00000351666               0    core     -1
#> 31   ensembl_homo_sapiens_37 ENST00000309486               0    core     -1
#>    assembly_name seq_region_name species         source      end
#> 1         GRCh37              17   human ensembl_havana 41277387
#> 2         GRCh37              17   human ensembl_havana 41277468
#> 3         GRCh37              17   human         havana 41277346
#> 4         GRCh37              17   human         havana 41277346
#> 5         GRCh37              17   human         havana 41277346
#> 6         GRCh37              17   human ensembl_havana 41277419
#> 7         GRCh37              17   human ensembl_havana 41277500
#> 8         GRCh37              17   human         havana 41277305
#> 9         GRCh37              17   human         havana 41277373
#> 10        GRCh37              17   human         havana 41256933
#> 11        GRCh37              17   human         havana 41277376
#> 12        GRCh37              17   human         havana 41277419
#> 13        GRCh37              17   human         havana 41223083
#> 14        GRCh37              17   human         havana 41256933
#> 15        GRCh37              17   human         havana 41243841
#> 16        GRCh37              17   human         havana 41277332
#> 17        GRCh37              17   human         havana 41247883
#> 18        GRCh37              17   human         havana 41277340
#> 19        GRCh37              17   human         havana 41277381
#> 20        GRCh37              17   human         havana 41277317
#> 21        GRCh37              17   human         havana 41277317
#> 22        GRCh37              17   human         havana 41277467
#> 23        GRCh37              17   human         havana 41256908
#> 24        GRCh37              17   human         havana 41277370
#> 25        GRCh37              17   human         havana 41277387
#> 26        GRCh37              17   human         havana 41277338
#> 27        GRCh37              17   human        ensembl 41277500
#> 28        GRCh37              17   human        ensembl 41277500
#> 29        GRCh37              17   human        ensembl 41277468
#> 30        GRCh37              17   human        ensembl 41276132
#> 31        GRCh37              17   human        ensembl 41277468
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $version
#> [1] 10
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $Transcript
#>   db_type      end          Parent seq_region_name      species         Exon
#> 1    core 32973347 ENSG00000139618              13 homo_sapiens c(1, 1, ....
#> 2    core 32907428 ENSG00000139618              13 homo_sapiens c(1, 1, ....
#> 3    core 32953632 ENSG00000139618              13 homo_sapiens c(329451....
#> 4    core 32972409 ENSG00000139618              13 homo_sapiens c("ENSE0....
#> 5    core 32972585 ENSG00000139618              13 homo_sapiens c(329709....
#> 6    core 32973805 ENSG00000139618              13 homo_sapiens c("GRCh3....
#>   is_canonical strand                logic_name                 biotype
#> 1            0      1 ensembl_havana_transcript          protein_coding
#> 2            0      1    havana_homo_sapiens_37          protein_coding
#> 3            0      1    havana_homo_sapiens_37 nonsense_mediated_decay
#> 4            0      1    havana_homo_sapiens_37 nonsense_mediated_decay
#> 5            0      1    havana_homo_sapiens_37         retained_intron
#> 6            1      1   ensembl_homo_sapiens_37          protein_coding
#>   object_type version assembly_name              id gencode_primary
#> 1  Transcript       3        GRCh37 ENST00000380152               0
#> 2  Transcript       2        GRCh37 ENST00000530893               0
#> 3  Transcript       1        GRCh37 ENST00000528762               0
#> 4  Transcript       1        GRCh37 ENST00000470094               0
#> 5  Transcript       1        GRCh37 ENST00000533776               0
#> 6  Transcript       1        GRCh37 ENST00000544455               0
#>   display_name Translation.version Translation.db_type Translation.object_type
#> 1    BRCA2-001                   3                core             Translation
#> 2    BRCA2-003                   2                core             Translation
#> 3    BRCA2-005                   1                core             Translation
#> 4    BRCA2-002                   1                core             Translation
#> 5    BRCA2-006                  NA                <NA>                    <NA>
#> 6    BRCA2-201                   1                core             Translation
#>   Translation.Parent Translation.end  Translation.id Translation.start
#> 1    ENST00000380152        32972907 ENSP00000369497          32890598
#> 2    ENST00000530893        32907428 ENSP00000435699          32899266
#> 3    ENST00000528762        32950807 ENSP00000433168          32945108
#> 4    ENST00000470094        32970229 ENSP00000434898          32953977
#> 5               <NA>              NA            <NA>                NA
#> 6    ENST00000544455        32972907 ENSP00000439902          32890598
#>   Translation.species Translation.length    start length         source
#> 1        homo_sapiens               3418 32889611  10930 ensembl_havana
#> 2        homo_sapiens                481 32889642   2011         havana
#> 3        homo_sapiens                 64 32945108    495         havana
#> 4        homo_sapiens                186 32953977    842         havana
#> 5                <NA>                 NA 32970946    523         havana
#> 6        homo_sapiens               3418 32889617  10984        ensembl
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $start
#> [1] 32889611
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $db_type
#> [1] "core"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $end
#> [1] 32973805
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $strand
#> [1] 1
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $version
#> [1] 15
#> 
#> $Transcript
#>                 id assembly_name display_name    start length gencode_primary
#> 1  ENST00000357654        GRCh37    BRCA1-001 41196312   7094               0
#> 2  ENST00000468300        GRCh37    BRCA1-007 41196822   3273               0
#> 3  ENST00000586385        GRCh37    BRCA1-023 41197580    781               0
#> 4  ENST00000591534        GRCh37    BRCA1-024 41197580   1282               0
#> 5  ENST00000591849        GRCh37    BRCA1-025 41197580    563               0
#> 6  ENST00000493795        GRCh37    BRCA1-006 41197646   5732               0
#> 7  ENST00000471181        GRCh37    BRCA1-005 41197646   5936               0
#> 8  ENST00000461221        GRCh37    BRCA1-010 41197695   5693               0
#> 9  ENST00000491747        GRCh37    BRCA1-014 41197695   2379               0
#> 10 ENST00000484087        GRCh37    BRCA1-015 41215361   1495               0
#> 11 ENST00000478531        GRCh37    BRCA1-009 41215361   1972               0
#> 12 ENST00000493919        GRCh37    BRCA1-008 41215377   1948               0
#> 13 ENST00000472490        GRCh37    BRCA1-021 41219291    561               0
#> 14 ENST00000487825        GRCh37    BRCA1-019 41228505    800               0
#> 15 ENST00000461574        GRCh37    BRCA1-022 41228554    726               0
#> 16 ENST00000467274        GRCh37    BRCA1-012 41243115   4497               0
#> 17 ENST00000412061        GRCh37    BRCA1-026 41245587   1312               0
#> 18 ENST00000470026        GRCh37    BRCA1-011 41245601   2108               0
#> 19 ENST00000477152        GRCh37    BRCA1-004 41245603   1980               0
#> 20 ENST00000492859        GRCh37    BRCA1-002 41246129   1584               0
#> 21 ENST00000497488        GRCh37    BRCA1-003 41246129    779               0
#> 22 ENST00000494123        GRCh37    BRCA1-013 41246129   1612               0
#> 23 ENST00000473961        GRCh37    BRCA1-018 41246187    958               0
#> 24 ENST00000476777        GRCh37    BRCA1-017 41247863    769               0
#> 25 ENST00000461798        GRCh37    BRCA1-020 41251848    582               0
#> 26 ENST00000489037        GRCh37    BRCA1-016 41256206    455               0
#> 27 ENST00000354071        GRCh37    BRCA1-205 41196313   6411               0
#> 28 ENST00000352993        GRCh37    BRCA1-204 41196313   3780               0
#> 29 ENST00000346315        GRCh37    BRCA1-202 41196313   6451               0
#> 30 ENST00000351666        GRCh37    BRCA1-203 41196313   3444               0
#> 31 ENST00000309486        GRCh37    BRCA1-201 41196313   7114               0
#>    version db_type         Exon object_type Translation.end  Translation.id
#> 1        3    core c(412773....  Transcript        41276113 ENSP00000350283
#> 2        1    core c("Exon"....  Transcript        41276113 ENSP00000417148
#> 3        1    core c("core"....  Transcript        41277202 ENSP00000465818
#> 4        1    core c(412772....  Transcript        41226495 ENSP00000467329
#> 5        1    core c(412773....  Transcript        41202109 ENSP00000465347
#> 6        1    core c("Exon"....  Transcript        41258543 ENSP00000418775
#> 7        2    core c(412772....  Transcript        41276113 ENSP00000418960
#> 8        1    core c("Exon"....  Transcript        41276113 ENSP00000418548
#> 9        2    core c(412772....  Transcript        41276113 ENSP00000420705
#> 10       1    core c("Exon"....  Transcript        41256933 ENSP00000419481
#> 11       1    core c("Exon"....  Transcript        41276113 ENSP00000420412
#> 12       1    core c("GRCh3....  Transcript        41258543 ENSP00000418819
#> 13       1    core c("ENSE0....  Transcript              NA            <NA>
#> 14       1    core c("human....  Transcript        41256933 ENSP00000418212
#> 15       1    core c("Exon"....  Transcript        41243841 ENSP00000417241
#> 16       1    core c("Exon"....  Transcript              NA            <NA>
#> 17       3    core c("Exon"....  Transcript        41247883 ENSP00000397145
#> 18       1    core c("core"....  Transcript        41276113 ENSP00000419274
#> 19       1    core c("core"....  Transcript        41276113 ENSP00000419988
#> 20       1    core c("ENSE0....  Transcript        41276113 ENSP00000420253
#> 21       1    core c("Exon"....  Transcript        41246659 ENSP00000418986
#> 22       1    core c("Exon"....  Transcript        41276113 ENSP00000419103
#> 23       1    core c("Exon"....  Transcript        41256908 ENSP00000420201
#> 24       1    core c("core"....  Transcript        41276113 ENSP00000417554
#> 25       1    core c("Exon"....  Transcript        41276113 ENSP00000417988
#> 26       1    core c("ENSE0....  Transcript        41276113 ENSP00000420781
#> 27       3    core c("Exon"....  Transcript        41276113 ENSP00000326002
#> 28       3    core c("core"....  Transcript        41276113 ENSP00000312236
#> 29       3    core c("Exon"....  Transcript        41276113 ENSP00000246907
#> 30       3    core c("core"....  Transcript        41276113 ENSP00000338007
#> 31       4    core c("core"....  Transcript        41246659 ENSP00000310938
#>    Translation.species Translation.version Translation.length
#> 1                human                   3               1863
#> 2                human                   1                699
#> 3                human                   1                173
#> 4                human                   1                354
#> 5                human                   1                 96
#> 6                human                   1               1816
#> 7                human                   2               1884
#> 8                human                   1                 63
#> 9                human                   2                759
#> 10               human                   1                498
#> 11               human                   1                623
#> 12               human                   1                572
#> 13                <NA>                  NA                 NA
#> 14               human                   1                266
#> 15               human                   1                242
#> 16                <NA>                  NA                 NA
#> 17               human                   3                437
#> 18               human                   1                649
#> 19               human                   1                622
#> 20               human                   1                 59
#> 21               human                   1                177
#> 22               human                   1                473
#> 23               human                   1                319
#> 24               human                   1                222
#> 25               human                   1                 63
#> 26               human                   1                 98
#> 27               human                   6               1598
#> 28               human                   5                721
#> 29               human                   4               1624
#> 30               human                   3                680
#> 31               human                   4               1567
#>    Translation.Parent Translation.start Translation.db_type
#> 1     ENST00000357654          41197695                core
#> 2     ENST00000468300          41197801                core
#> 3     ENST00000586385          41197695                core
#> 4     ENST00000591534          41197695                core
#> 5     ENST00000591849          41197695                core
#> 6     ENST00000493795          41197695                core
#> 7     ENST00000471181          41197695                core
#> 8     ENST00000461221          41256972                core
#> 9     ENST00000491747          41197695                core
#> 10    ENST00000484087          41215361                core
#> 11    ENST00000478531          41215361                core
#> 12    ENST00000493919          41215377                core
#> 13               <NA>                NA                <NA>
#> 14    ENST00000487825          41228505                core
#> 15    ENST00000461574          41228554                core
#> 16               <NA>                NA                <NA>
#> 17    ENST00000412061          41245587                core
#> 18    ENST00000470026          41245601                core
#> 19    ENST00000477152          41245603                core
#> 20    ENST00000492859          41262552                core
#> 21    ENST00000497488          41246129                core
#> 22    ENST00000494123          41246129                core
#> 23    ENST00000473961          41246187                core
#> 24    ENST00000476777          41247863                core
#> 25    ENST00000461798          41256972                core
#> 26    ENST00000489037          41256206                core
#> 27    ENST00000354071          41197695                core
#> 28    ENST00000352993          41197695                core
#> 29    ENST00000346315          41197695                core
#> 30    ENST00000351666          41197695                core
#> 31    ENST00000309486          41197695                core
#>    Translation.object_type strand                logic_name      end
#> 1              Translation     -1 ensembl_havana_transcript 41277387
#> 2              Translation     -1 ensembl_havana_transcript 41277468
#> 3              Translation     -1    havana_homo_sapiens_37 41277346
#> 4              Translation     -1    havana_homo_sapiens_37 41277346
#> 5              Translation     -1    havana_homo_sapiens_37 41277346
#> 6              Translation     -1 ensembl_havana_transcript 41277419
#> 7              Translation     -1 ensembl_havana_transcript 41277500
#> 8              Translation     -1    havana_homo_sapiens_37 41277305
#> 9              Translation     -1    havana_homo_sapiens_37 41277373
#> 10             Translation     -1    havana_homo_sapiens_37 41256933
#> 11             Translation     -1    havana_homo_sapiens_37 41277376
#> 12             Translation     -1    havana_homo_sapiens_37 41277419
#> 13                    <NA>     -1    havana_homo_sapiens_37 41223083
#> 14             Translation     -1    havana_homo_sapiens_37 41256933
#> 15             Translation     -1    havana_homo_sapiens_37 41243841
#> 16                    <NA>     -1    havana_homo_sapiens_37 41277332
#> 17             Translation     -1    havana_homo_sapiens_37 41247883
#> 18             Translation     -1    havana_homo_sapiens_37 41277340
#> 19             Translation     -1    havana_homo_sapiens_37 41277381
#> 20             Translation     -1    havana_homo_sapiens_37 41277317
#> 21             Translation     -1    havana_homo_sapiens_37 41277317
#> 22             Translation     -1    havana_homo_sapiens_37 41277467
#> 23             Translation     -1    havana_homo_sapiens_37 41256908
#> 24             Translation     -1    havana_homo_sapiens_37 41277370
#> 25             Translation     -1    havana_homo_sapiens_37 41277387
#> 26             Translation     -1    havana_homo_sapiens_37 41277338
#> 27             Translation     -1   ensembl_homo_sapiens_37 41277500
#> 28             Translation     -1   ensembl_homo_sapiens_37 41277500
#> 29             Translation     -1   ensembl_homo_sapiens_37 41277468
#> 30             Translation     -1   ensembl_homo_sapiens_37 41276132
#> 31             Translation     -1   ensembl_homo_sapiens_37 41277468
#>             Parent seq_region_name                 biotype species is_canonical
#> 1  ENSG00000012048              17          protein_coding   human            0
#> 2  ENSG00000012048              17          protein_coding   human            0
#> 3  ENSG00000012048              17          protein_coding   human            0
#> 4  ENSG00000012048              17          protein_coding   human            0
#> 5  ENSG00000012048              17          protein_coding   human            0
#> 6  ENSG00000012048              17          protein_coding   human            0
#> 7  ENSG00000012048              17          protein_coding   human            1
#> 8  ENSG00000012048              17 nonsense_mediated_decay   human            0
#> 9  ENSG00000012048              17          protein_coding   human            0
#> 10 ENSG00000012048              17          protein_coding   human            0
#> 11 ENSG00000012048              17          protein_coding   human            0
#> 12 ENSG00000012048              17          protein_coding   human            0
#> 13 ENSG00000012048              17         retained_intron   human            0
#> 14 ENSG00000012048              17          protein_coding   human            0
#> 15 ENSG00000012048              17          protein_coding   human            0
#> 16 ENSG00000012048              17         retained_intron   human            0
#> 17 ENSG00000012048              17          non_stop_decay   human            0
#> 18 ENSG00000012048              17          protein_coding   human            0
#> 19 ENSG00000012048              17          protein_coding   human            0
#> 20 ENSG00000012048              17 nonsense_mediated_decay   human            0
#> 21 ENSG00000012048              17          protein_coding   human            0
#> 22 ENSG00000012048              17          protein_coding   human            0
#> 23 ENSG00000012048              17          protein_coding   human            0
#> 24 ENSG00000012048              17          protein_coding   human            0
#> 25 ENSG00000012048              17 nonsense_mediated_decay   human            0
#> 26 ENSG00000012048              17          protein_coding   human            0
#> 27 ENSG00000012048              17          protein_coding   human            0
#> 28 ENSG00000012048              17          protein_coding   human            0
#> 29 ENSG00000012048              17          protein_coding   human            0
#> 30 ENSG00000012048              17          protein_coding   human            0
#> 31 ENSG00000012048              17          protein_coding   human            0
#>            source
#> 1  ensembl_havana
#> 2  ensembl_havana
#> 3          havana
#> 4          havana
#> 5          havana
#> 6  ensembl_havana
#> 7  ensembl_havana
#> 8          havana
#> 9          havana
#> 10         havana
#> 11         havana
#> 12         havana
#> 13         havana
#> 14         havana
#> 15         havana
#> 16         havana
#> 17         havana
#> 18         havana
#> 19         havana
#> 20         havana
#> 21         havana
#> 22         havana
#> 23         havana
#> 24         havana
#> 25         havana
#> 26         havana
#> 27        ensembl
#> 28        ensembl
#> 29        ensembl
#> 30        ensembl
#> 31        ensembl
#> 
#> $start
#> [1] 41196312
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $species
#> [1] "human"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $end
#> [1] 41277500
#> 
#> $strand
#> [1] -1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
```
