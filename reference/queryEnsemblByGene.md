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
#> $source
#> [1] "ensembl_havana"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $start
#> [1] 41196312
#> 
#> $Transcript
#>    gencode_primary              id assembly_name version
#> 1                0 ENST00000357654        GRCh37       3
#> 2                0 ENST00000468300        GRCh37       1
#> 3                0 ENST00000586385        GRCh37       1
#> 4                0 ENST00000591534        GRCh37       1
#> 5                0 ENST00000591849        GRCh37       1
#> 6                0 ENST00000493795        GRCh37       1
#> 7                0 ENST00000471181        GRCh37       2
#> 8                0 ENST00000461221        GRCh37       1
#> 9                0 ENST00000491747        GRCh37       2
#> 10               0 ENST00000484087        GRCh37       1
#> 11               0 ENST00000478531        GRCh37       1
#> 12               0 ENST00000493919        GRCh37       1
#> 13               0 ENST00000472490        GRCh37       1
#> 14               0 ENST00000487825        GRCh37       1
#> 15               0 ENST00000461574        GRCh37       1
#> 16               0 ENST00000467274        GRCh37       1
#> 17               0 ENST00000412061        GRCh37       3
#> 18               0 ENST00000470026        GRCh37       1
#> 19               0 ENST00000477152        GRCh37       1
#> 20               0 ENST00000492859        GRCh37       1
#> 21               0 ENST00000497488        GRCh37       1
#> 22               0 ENST00000494123        GRCh37       1
#> 23               0 ENST00000473961        GRCh37       1
#> 24               0 ENST00000476777        GRCh37       1
#> 25               0 ENST00000461798        GRCh37       1
#> 26               0 ENST00000489037        GRCh37       1
#> 27               0 ENST00000354071        GRCh37       3
#> 28               0 ENST00000352993        GRCh37       3
#> 29               0 ENST00000346315        GRCh37       3
#> 30               0 ENST00000351666        GRCh37       3
#> 31               0 ENST00000309486        GRCh37       4
#>                   logic_name                 biotype object_type length
#> 1  ensembl_havana_transcript          protein_coding  Transcript   7094
#> 2  ensembl_havana_transcript          protein_coding  Transcript   3273
#> 3     havana_homo_sapiens_37          protein_coding  Transcript    781
#> 4     havana_homo_sapiens_37          protein_coding  Transcript   1282
#> 5     havana_homo_sapiens_37          protein_coding  Transcript    563
#> 6  ensembl_havana_transcript          protein_coding  Transcript   5732
#> 7  ensembl_havana_transcript          protein_coding  Transcript   5936
#> 8     havana_homo_sapiens_37 nonsense_mediated_decay  Transcript   5693
#> 9     havana_homo_sapiens_37          protein_coding  Transcript   2379
#> 10    havana_homo_sapiens_37          protein_coding  Transcript   1495
#> 11    havana_homo_sapiens_37          protein_coding  Transcript   1972
#> 12    havana_homo_sapiens_37          protein_coding  Transcript   1948
#> 13    havana_homo_sapiens_37         retained_intron  Transcript    561
#> 14    havana_homo_sapiens_37          protein_coding  Transcript    800
#> 15    havana_homo_sapiens_37          protein_coding  Transcript    726
#> 16    havana_homo_sapiens_37         retained_intron  Transcript   4497
#> 17    havana_homo_sapiens_37          non_stop_decay  Transcript   1312
#> 18    havana_homo_sapiens_37          protein_coding  Transcript   2108
#> 19    havana_homo_sapiens_37          protein_coding  Transcript   1980
#> 20    havana_homo_sapiens_37 nonsense_mediated_decay  Transcript   1584
#> 21    havana_homo_sapiens_37          protein_coding  Transcript    779
#> 22    havana_homo_sapiens_37          protein_coding  Transcript   1612
#> 23    havana_homo_sapiens_37          protein_coding  Transcript    958
#> 24    havana_homo_sapiens_37          protein_coding  Transcript    769
#> 25    havana_homo_sapiens_37 nonsense_mediated_decay  Transcript    582
#> 26    havana_homo_sapiens_37          protein_coding  Transcript    455
#> 27   ensembl_homo_sapiens_37          protein_coding  Transcript   6411
#> 28   ensembl_homo_sapiens_37          protein_coding  Transcript   3780
#> 29   ensembl_homo_sapiens_37          protein_coding  Transcript   6451
#> 30   ensembl_homo_sapiens_37          protein_coding  Transcript   3444
#> 31   ensembl_homo_sapiens_37          protein_coding  Transcript   7114
#>            source    start display_name Translation.version Translation.db_type
#> 1  ensembl_havana 41196312    BRCA1-001                   3                core
#> 2  ensembl_havana 41196822    BRCA1-007                   1                core
#> 3          havana 41197580    BRCA1-023                   1                core
#> 4          havana 41197580    BRCA1-024                   1                core
#> 5          havana 41197580    BRCA1-025                   1                core
#> 6  ensembl_havana 41197646    BRCA1-006                   1                core
#> 7  ensembl_havana 41197646    BRCA1-005                   2                core
#> 8          havana 41197695    BRCA1-010                   1                core
#> 9          havana 41197695    BRCA1-014                   2                core
#> 10         havana 41215361    BRCA1-015                   1                core
#> 11         havana 41215361    BRCA1-009                   1                core
#> 12         havana 41215377    BRCA1-008                   1                core
#> 13         havana 41219291    BRCA1-021                  NA                <NA>
#> 14         havana 41228505    BRCA1-019                   1                core
#> 15         havana 41228554    BRCA1-022                   1                core
#> 16         havana 41243115    BRCA1-012                  NA                <NA>
#> 17         havana 41245587    BRCA1-026                   3                core
#> 18         havana 41245601    BRCA1-011                   1                core
#> 19         havana 41245603    BRCA1-004                   1                core
#> 20         havana 41246129    BRCA1-002                   1                core
#> 21         havana 41246129    BRCA1-003                   1                core
#> 22         havana 41246129    BRCA1-013                   1                core
#> 23         havana 41246187    BRCA1-018                   1                core
#> 24         havana 41247863    BRCA1-017                   1                core
#> 25         havana 41251848    BRCA1-020                   1                core
#> 26         havana 41256206    BRCA1-016                   1                core
#> 27        ensembl 41196313    BRCA1-205                   6                core
#> 28        ensembl 41196313    BRCA1-204                   5                core
#> 29        ensembl 41196313    BRCA1-202                   4                core
#> 30        ensembl 41196313    BRCA1-203                   3                core
#> 31        ensembl 41196313    BRCA1-201                   4                core
#>    Translation.object_type Translation.Parent Translation.end  Translation.id
#> 1              Translation    ENST00000357654        41276113 ENSP00000350283
#> 2              Translation    ENST00000468300        41276113 ENSP00000417148
#> 3              Translation    ENST00000586385        41277202 ENSP00000465818
#> 4              Translation    ENST00000591534        41226495 ENSP00000467329
#> 5              Translation    ENST00000591849        41202109 ENSP00000465347
#> 6              Translation    ENST00000493795        41258543 ENSP00000418775
#> 7              Translation    ENST00000471181        41276113 ENSP00000418960
#> 8              Translation    ENST00000461221        41276113 ENSP00000418548
#> 9              Translation    ENST00000491747        41276113 ENSP00000420705
#> 10             Translation    ENST00000484087        41256933 ENSP00000419481
#> 11             Translation    ENST00000478531        41276113 ENSP00000420412
#> 12             Translation    ENST00000493919        41258543 ENSP00000418819
#> 13                    <NA>               <NA>              NA            <NA>
#> 14             Translation    ENST00000487825        41256933 ENSP00000418212
#> 15             Translation    ENST00000461574        41243841 ENSP00000417241
#> 16                    <NA>               <NA>              NA            <NA>
#> 17             Translation    ENST00000412061        41247883 ENSP00000397145
#> 18             Translation    ENST00000470026        41276113 ENSP00000419274
#> 19             Translation    ENST00000477152        41276113 ENSP00000419988
#> 20             Translation    ENST00000492859        41276113 ENSP00000420253
#> 21             Translation    ENST00000497488        41246659 ENSP00000418986
#> 22             Translation    ENST00000494123        41276113 ENSP00000419103
#> 23             Translation    ENST00000473961        41256908 ENSP00000420201
#> 24             Translation    ENST00000476777        41276113 ENSP00000417554
#> 25             Translation    ENST00000461798        41276113 ENSP00000417988
#> 26             Translation    ENST00000489037        41276113 ENSP00000420781
#> 27             Translation    ENST00000354071        41276113 ENSP00000326002
#> 28             Translation    ENST00000352993        41276113 ENSP00000312236
#> 29             Translation    ENST00000346315        41276113 ENSP00000246907
#> 30             Translation    ENST00000351666        41276113 ENSP00000338007
#> 31             Translation    ENST00000309486        41246659 ENSP00000310938
#>    Translation.start Translation.species Translation.length          Parent
#> 1           41197695               human               1863 ENSG00000012048
#> 2           41197801               human                699 ENSG00000012048
#> 3           41197695               human                173 ENSG00000012048
#> 4           41197695               human                354 ENSG00000012048
#> 5           41197695               human                 96 ENSG00000012048
#> 6           41197695               human               1816 ENSG00000012048
#> 7           41197695               human               1884 ENSG00000012048
#> 8           41256972               human                 63 ENSG00000012048
#> 9           41197695               human                759 ENSG00000012048
#> 10          41215361               human                498 ENSG00000012048
#> 11          41215361               human                623 ENSG00000012048
#> 12          41215377               human                572 ENSG00000012048
#> 13                NA                <NA>                 NA ENSG00000012048
#> 14          41228505               human                266 ENSG00000012048
#> 15          41228554               human                242 ENSG00000012048
#> 16                NA                <NA>                 NA ENSG00000012048
#> 17          41245587               human                437 ENSG00000012048
#> 18          41245601               human                649 ENSG00000012048
#> 19          41245603               human                622 ENSG00000012048
#> 20          41262552               human                 59 ENSG00000012048
#> 21          41246129               human                177 ENSG00000012048
#> 22          41246129               human                473 ENSG00000012048
#> 23          41246187               human                319 ENSG00000012048
#> 24          41247863               human                222 ENSG00000012048
#> 25          41256972               human                 63 ENSG00000012048
#> 26          41256206               human                 98 ENSG00000012048
#> 27          41197695               human               1598 ENSG00000012048
#> 28          41197695               human                721 ENSG00000012048
#> 29          41197695               human               1624 ENSG00000012048
#> 30          41197695               human                680 ENSG00000012048
#> 31          41197695               human               1567 ENSG00000012048
#>    seq_region_name      end db_type strand is_canonical         Exon species
#> 1               17 41277387    core     -1            0 c("human....   human
#> 2               17 41277468    core     -1            0 c("human....   human
#> 3               17 41277346    core     -1            0 c("17", ....   human
#> 4               17 41277346    core     -1            0 c("17", ....   human
#> 5               17 41277346    core     -1            0 c(-1, -1....   human
#> 6               17 41277419    core     -1            0 c(1, 1, ....   human
#> 7               17 41277500    core     -1            1 c(412775....   human
#> 8               17 41277305    core     -1            0 c(412771....   human
#> 9               17 41277373    core     -1            0 c(-1, -1....   human
#> 10              17 41256933    core     -1            0 c(-1, -1....   human
#> 11              17 41277376    core     -1            0 c("17", ....   human
#> 12              17 41277419    core     -1            0 c("core"....   human
#> 13              17 41223083    core     -1            0 c(1, 1),....   human
#> 14              17 41256933    core     -1            0 c(-1, -1....   human
#> 15              17 41243841    core     -1            0 c("ENSE0....   human
#> 16              17 41277332    core     -1            0 c(-1, -1....   human
#> 17              17 41247883    core     -1            0 c("17", ....   human
#> 18              17 41277340    core     -1            0 c("17", ....   human
#> 19              17 41277381    core     -1            0 c("human....   human
#> 20              17 41277317    core     -1            0 c(-1, -1....   human
#> 21              17 41277317    core     -1            0 c(1, 1),....   human
#> 22              17 41277467    core     -1            0 c("GRCh3....   human
#> 23              17 41256908    core     -1            0 c(412569....   human
#> 24              17 41277370    core     -1            0 c(-1, -1....   human
#> 25              17 41277387    core     -1            0 c(-1, -1....   human
#> 26              17 41277338    core     -1            0 c(-1, -1....   human
#> 27              17 41277500    core     -1            0 c(-1, -1....   human
#> 28              17 41277500    core     -1            0 c("core"....   human
#> 29              17 41277468    core     -1            0 c("core"....   human
#> 30              17 41276132    core     -1            0 c("core"....   human
#> 31              17 41277468    core     -1            0 c("17", ....   human
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
#> $strand
#> [1] -1
#> 
#> $species
#> [1] "human"
#> 
#> $end
#> [1] 41277500
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $db_type
#> [1] "core"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $source
#> [1] "ensembl_havana"
#> 
#> $strand
#> [1] 1
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
#> $Transcript
#>      start length display_name gencode_primary Translation.length
#> 1 32889611  10930    BRCA2-001               0               3418
#> 2 32889642   2011    BRCA2-003               0                481
#> 3 32945108    495    BRCA2-005               0                 64
#> 4 32953977    842    BRCA2-002               0                186
#> 5 32970946    523    BRCA2-006               0                 NA
#> 6 32889617  10984    BRCA2-201               0               3418
#>   Translation.start  Translation.id Translation.version Translation.db_type
#> 1          32890598 ENSP00000369497                   3                core
#> 2          32899266 ENSP00000435699                   2                core
#> 3          32945108 ENSP00000433168                   1                core
#> 4          32953977 ENSP00000434898                   1                core
#> 5                NA            <NA>                  NA                <NA>
#> 6          32890598 ENSP00000439902                   1                core
#>   Translation.end Translation.Parent Translation.object_type
#> 1        32972907    ENST00000380152             Translation
#> 2        32907428    ENST00000530893             Translation
#> 3        32950807    ENST00000528762             Translation
#> 4        32970229    ENST00000470094             Translation
#> 5              NA               <NA>                    <NA>
#> 6        32972907    ENST00000544455             Translation
#>   Translation.species      end version assembly_name db_type          Parent
#> 1        homo_sapiens 32973347       3        GRCh37    core ENSG00000139618
#> 2        homo_sapiens 32907428       2        GRCh37    core ENSG00000139618
#> 3        homo_sapiens 32953632       1        GRCh37    core ENSG00000139618
#> 4        homo_sapiens 32972409       1        GRCh37    core ENSG00000139618
#> 5                <NA> 32972585       1        GRCh37    core ENSG00000139618
#> 6        homo_sapiens 32973805       1        GRCh37    core ENSG00000139618
#>   seq_region_name      species strand         source                logic_name
#> 1              13 homo_sapiens      1 ensembl_havana ensembl_havana_transcript
#> 2              13 homo_sapiens      1         havana    havana_homo_sapiens_37
#> 3              13 homo_sapiens      1         havana    havana_homo_sapiens_37
#> 4              13 homo_sapiens      1         havana    havana_homo_sapiens_37
#> 5              13 homo_sapiens      1         havana    havana_homo_sapiens_37
#> 6              13 homo_sapiens      1        ensembl   ensembl_homo_sapiens_37
#>   is_canonical              id                 biotype         Exon object_type
#> 1            0 ENST00000380152          protein_coding c(328898....  Transcript
#> 2            0 ENST00000530893          protein_coding c("Exon"....  Transcript
#> 3            0 ENST00000528762 nonsense_mediated_decay c(329452....  Transcript
#> 4            0 ENST00000470094 nonsense_mediated_decay c(1, 1, ....  Transcript
#> 5            0 ENST00000533776         retained_intron c(1, 1),....  Transcript
#> 6            1 ENST00000544455          protein_coding c("13", ....  Transcript
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $start
#> [1] 32889611
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $end
#> [1] 32973805
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $version
#> [1] 10
#> 
#> $db_type
#> [1] "core"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $biotype
#> [1] "protein_coding"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $start
#> [1] 41196312
#> 
#> $end
#> [1] 41277500
#> 
#> $Transcript
#>    gencode_primary                logic_name         Exon db_type      end
#> 1                0 ensembl_havana_transcript c(412772....    core 41277387
#> 2                0 ensembl_havana_transcript c(-1, -1....    core 41277468
#> 3                0    havana_homo_sapiens_37 c(412773....    core 41277346
#> 4                0    havana_homo_sapiens_37 c("human....    core 41277346
#> 5                0    havana_homo_sapiens_37 c("human....    core 41277346
#> 6                0 ensembl_havana_transcript c("ENSE0....    core 41277419
#> 7                0 ensembl_havana_transcript c("ENSE0....    core 41277500
#> 8                0    havana_homo_sapiens_37 c(-1, -1....    core 41277305
#> 9                0    havana_homo_sapiens_37 c(1, 1, ....    core 41277373
#> 10               0    havana_homo_sapiens_37 c("Exon"....    core 41256933
#> 11               0    havana_homo_sapiens_37 c("Exon"....    core 41277376
#> 12               0    havana_homo_sapiens_37 c(-1, -1....    core 41277419
#> 13               0    havana_homo_sapiens_37 c("human....    core 41223083
#> 14               0    havana_homo_sapiens_37 c("core"....    core 41256933
#> 15               0    havana_homo_sapiens_37 c(-1, -1....    core 41243841
#> 16               0    havana_homo_sapiens_37 c(1, 1, ....    core 41277332
#> 17               0    havana_homo_sapiens_37 c("Exon"....    core 41247883
#> 18               0    havana_homo_sapiens_37 c(1, 1, ....    core 41277340
#> 19               0    havana_homo_sapiens_37 c("ENSE0....    core 41277381
#> 20               0    havana_homo_sapiens_37 c(412773....    core 41277317
#> 21               0    havana_homo_sapiens_37 c("human....    core 41277317
#> 22               0    havana_homo_sapiens_37 c(412774....    core 41277467
#> 23               0    havana_homo_sapiens_37 c(1, 1, ....    core 41256908
#> 24               0    havana_homo_sapiens_37 c("human....    core 41277370
#> 25               0    havana_homo_sapiens_37 c(1, 1, ....    core 41277387
#> 26               0    havana_homo_sapiens_37 c(412773....    core 41277338
#> 27               0   ensembl_homo_sapiens_37 c(1, 1, ....    core 41277500
#> 28               0   ensembl_homo_sapiens_37 c(1, 1, ....    core 41277500
#> 29               0   ensembl_homo_sapiens_37 c("17", ....    core 41277468
#> 30               0   ensembl_homo_sapiens_37 c("core"....    core 41276132
#> 31               0   ensembl_homo_sapiens_37 c(1, 1, ....    core 41277468
#>       start object_type length                 biotype is_canonical species
#> 1  41196312  Transcript   7094          protein_coding            0   human
#> 2  41196822  Transcript   3273          protein_coding            0   human
#> 3  41197580  Transcript    781          protein_coding            0   human
#> 4  41197580  Transcript   1282          protein_coding            0   human
#> 5  41197580  Transcript    563          protein_coding            0   human
#> 6  41197646  Transcript   5732          protein_coding            0   human
#> 7  41197646  Transcript   5936          protein_coding            1   human
#> 8  41197695  Transcript   5693 nonsense_mediated_decay            0   human
#> 9  41197695  Transcript   2379          protein_coding            0   human
#> 10 41215361  Transcript   1495          protein_coding            0   human
#> 11 41215361  Transcript   1972          protein_coding            0   human
#> 12 41215377  Transcript   1948          protein_coding            0   human
#> 13 41219291  Transcript    561         retained_intron            0   human
#> 14 41228505  Transcript    800          protein_coding            0   human
#> 15 41228554  Transcript    726          protein_coding            0   human
#> 16 41243115  Transcript   4497         retained_intron            0   human
#> 17 41245587  Transcript   1312          non_stop_decay            0   human
#> 18 41245601  Transcript   2108          protein_coding            0   human
#> 19 41245603  Transcript   1980          protein_coding            0   human
#> 20 41246129  Transcript   1584 nonsense_mediated_decay            0   human
#> 21 41246129  Transcript    779          protein_coding            0   human
#> 22 41246129  Transcript   1612          protein_coding            0   human
#> 23 41246187  Transcript    958          protein_coding            0   human
#> 24 41247863  Transcript    769          protein_coding            0   human
#> 25 41251848  Transcript    582 nonsense_mediated_decay            0   human
#> 26 41256206  Transcript    455          protein_coding            0   human
#> 27 41196313  Transcript   6411          protein_coding            0   human
#> 28 41196313  Transcript   3780          protein_coding            0   human
#> 29 41196313  Transcript   6451          protein_coding            0   human
#> 30 41196313  Transcript   3444          protein_coding            0   human
#> 31 41196313  Transcript   7114          protein_coding            0   human
#>    strand version display_name assembly_name         source          Parent
#> 1      -1       3    BRCA1-001        GRCh37 ensembl_havana ENSG00000012048
#> 2      -1       1    BRCA1-007        GRCh37 ensembl_havana ENSG00000012048
#> 3      -1       1    BRCA1-023        GRCh37         havana ENSG00000012048
#> 4      -1       1    BRCA1-024        GRCh37         havana ENSG00000012048
#> 5      -1       1    BRCA1-025        GRCh37         havana ENSG00000012048
#> 6      -1       1    BRCA1-006        GRCh37 ensembl_havana ENSG00000012048
#> 7      -1       2    BRCA1-005        GRCh37 ensembl_havana ENSG00000012048
#> 8      -1       1    BRCA1-010        GRCh37         havana ENSG00000012048
#> 9      -1       2    BRCA1-014        GRCh37         havana ENSG00000012048
#> 10     -1       1    BRCA1-015        GRCh37         havana ENSG00000012048
#> 11     -1       1    BRCA1-009        GRCh37         havana ENSG00000012048
#> 12     -1       1    BRCA1-008        GRCh37         havana ENSG00000012048
#> 13     -1       1    BRCA1-021        GRCh37         havana ENSG00000012048
#> 14     -1       1    BRCA1-019        GRCh37         havana ENSG00000012048
#> 15     -1       1    BRCA1-022        GRCh37         havana ENSG00000012048
#> 16     -1       1    BRCA1-012        GRCh37         havana ENSG00000012048
#> 17     -1       3    BRCA1-026        GRCh37         havana ENSG00000012048
#> 18     -1       1    BRCA1-011        GRCh37         havana ENSG00000012048
#> 19     -1       1    BRCA1-004        GRCh37         havana ENSG00000012048
#> 20     -1       1    BRCA1-002        GRCh37         havana ENSG00000012048
#> 21     -1       1    BRCA1-003        GRCh37         havana ENSG00000012048
#> 22     -1       1    BRCA1-013        GRCh37         havana ENSG00000012048
#> 23     -1       1    BRCA1-018        GRCh37         havana ENSG00000012048
#> 24     -1       1    BRCA1-017        GRCh37         havana ENSG00000012048
#> 25     -1       1    BRCA1-020        GRCh37         havana ENSG00000012048
#> 26     -1       1    BRCA1-016        GRCh37         havana ENSG00000012048
#> 27     -1       3    BRCA1-205        GRCh37        ensembl ENSG00000012048
#> 28     -1       3    BRCA1-204        GRCh37        ensembl ENSG00000012048
#> 29     -1       3    BRCA1-202        GRCh37        ensembl ENSG00000012048
#> 30     -1       3    BRCA1-203        GRCh37        ensembl ENSG00000012048
#> 31     -1       4    BRCA1-201        GRCh37        ensembl ENSG00000012048
#>    seq_region_name Translation.Parent Translation.object_type  Translation.id
#> 1               17    ENST00000357654             Translation ENSP00000350283
#> 2               17    ENST00000468300             Translation ENSP00000417148
#> 3               17    ENST00000586385             Translation ENSP00000465818
#> 4               17    ENST00000591534             Translation ENSP00000467329
#> 5               17    ENST00000591849             Translation ENSP00000465347
#> 6               17    ENST00000493795             Translation ENSP00000418775
#> 7               17    ENST00000471181             Translation ENSP00000418960
#> 8               17    ENST00000461221             Translation ENSP00000418548
#> 9               17    ENST00000491747             Translation ENSP00000420705
#> 10              17    ENST00000484087             Translation ENSP00000419481
#> 11              17    ENST00000478531             Translation ENSP00000420412
#> 12              17    ENST00000493919             Translation ENSP00000418819
#> 13              17               <NA>                    <NA>            <NA>
#> 14              17    ENST00000487825             Translation ENSP00000418212
#> 15              17    ENST00000461574             Translation ENSP00000417241
#> 16              17               <NA>                    <NA>            <NA>
#> 17              17    ENST00000412061             Translation ENSP00000397145
#> 18              17    ENST00000470026             Translation ENSP00000419274
#> 19              17    ENST00000477152             Translation ENSP00000419988
#> 20              17    ENST00000492859             Translation ENSP00000420253
#> 21              17    ENST00000497488             Translation ENSP00000418986
#> 22              17    ENST00000494123             Translation ENSP00000419103
#> 23              17    ENST00000473961             Translation ENSP00000420201
#> 24              17    ENST00000476777             Translation ENSP00000417554
#> 25              17    ENST00000461798             Translation ENSP00000417988
#> 26              17    ENST00000489037             Translation ENSP00000420781
#> 27              17    ENST00000354071             Translation ENSP00000326002
#> 28              17    ENST00000352993             Translation ENSP00000312236
#> 29              17    ENST00000346315             Translation ENSP00000246907
#> 30              17    ENST00000351666             Translation ENSP00000338007
#> 31              17    ENST00000309486             Translation ENSP00000310938
#>    Translation.length Translation.end Translation.db_type Translation.start
#> 1                1863        41276113                core          41197695
#> 2                 699        41276113                core          41197801
#> 3                 173        41277202                core          41197695
#> 4                 354        41226495                core          41197695
#> 5                  96        41202109                core          41197695
#> 6                1816        41258543                core          41197695
#> 7                1884        41276113                core          41197695
#> 8                  63        41276113                core          41256972
#> 9                 759        41276113                core          41197695
#> 10                498        41256933                core          41215361
#> 11                623        41276113                core          41215361
#> 12                572        41258543                core          41215377
#> 13                 NA              NA                <NA>                NA
#> 14                266        41256933                core          41228505
#> 15                242        41243841                core          41228554
#> 16                 NA              NA                <NA>                NA
#> 17                437        41247883                core          41245587
#> 18                649        41276113                core          41245601
#> 19                622        41276113                core          41245603
#> 20                 59        41276113                core          41262552
#> 21                177        41246659                core          41246129
#> 22                473        41276113                core          41246129
#> 23                319        41256908                core          41246187
#> 24                222        41276113                core          41247863
#> 25                 63        41276113                core          41256972
#> 26                 98        41276113                core          41256206
#> 27               1598        41276113                core          41197695
#> 28                721        41276113                core          41197695
#> 29               1624        41276113                core          41197695
#> 30                680        41276113                core          41197695
#> 31               1567        41246659                core          41197695
#>    Translation.species Translation.version              id
#> 1                human                   3 ENST00000357654
#> 2                human                   1 ENST00000468300
#> 3                human                   1 ENST00000586385
#> 4                human                   1 ENST00000591534
#> 5                human                   1 ENST00000591849
#> 6                human                   1 ENST00000493795
#> 7                human                   2 ENST00000471181
#> 8                human                   1 ENST00000461221
#> 9                human                   2 ENST00000491747
#> 10               human                   1 ENST00000484087
#> 11               human                   1 ENST00000478531
#> 12               human                   1 ENST00000493919
#> 13                <NA>                  NA ENST00000472490
#> 14               human                   1 ENST00000487825
#> 15               human                   1 ENST00000461574
#> 16                <NA>                  NA ENST00000467274
#> 17               human                   3 ENST00000412061
#> 18               human                   1 ENST00000470026
#> 19               human                   1 ENST00000477152
#> 20               human                   1 ENST00000492859
#> 21               human                   1 ENST00000497488
#> 22               human                   1 ENST00000494123
#> 23               human                   1 ENST00000473961
#> 24               human                   1 ENST00000476777
#> 25               human                   1 ENST00000461798
#> 26               human                   1 ENST00000489037
#> 27               human                   6 ENST00000354071
#> 28               human                   5 ENST00000352993
#> 29               human                   4 ENST00000346315
#> 30               human                   3 ENST00000351666
#> 31               human                   4 ENST00000309486
#> 
#> $db_type
#> [1] "core"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $version
#> [1] 15
#> 
#> $strand
#> [1] -1
#> 
#> $species
#> [1] "human"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
```
