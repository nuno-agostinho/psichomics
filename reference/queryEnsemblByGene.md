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
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $Transcript
#>            source is_canonical          Parent seq_region_name species
#> 1  ensembl_havana            0 ENSG00000012048              17   human
#> 2  ensembl_havana            0 ENSG00000012048              17   human
#> 3          havana            0 ENSG00000012048              17   human
#> 4          havana            0 ENSG00000012048              17   human
#> 5          havana            0 ENSG00000012048              17   human
#> 6  ensembl_havana            0 ENSG00000012048              17   human
#> 7  ensembl_havana            1 ENSG00000012048              17   human
#> 8          havana            0 ENSG00000012048              17   human
#> 9          havana            0 ENSG00000012048              17   human
#> 10         havana            0 ENSG00000012048              17   human
#> 11         havana            0 ENSG00000012048              17   human
#> 12         havana            0 ENSG00000012048              17   human
#> 13         havana            0 ENSG00000012048              17   human
#> 14         havana            0 ENSG00000012048              17   human
#> 15         havana            0 ENSG00000012048              17   human
#> 16         havana            0 ENSG00000012048              17   human
#> 17         havana            0 ENSG00000012048              17   human
#> 18         havana            0 ENSG00000012048              17   human
#> 19         havana            0 ENSG00000012048              17   human
#> 20         havana            0 ENSG00000012048              17   human
#> 21         havana            0 ENSG00000012048              17   human
#> 22         havana            0 ENSG00000012048              17   human
#> 23         havana            0 ENSG00000012048              17   human
#> 24         havana            0 ENSG00000012048              17   human
#> 25         havana            0 ENSG00000012048              17   human
#> 26         havana            0 ENSG00000012048              17   human
#> 27        ensembl            0 ENSG00000012048              17   human
#> 28        ensembl            0 ENSG00000012048              17   human
#> 29        ensembl            0 ENSG00000012048              17   human
#> 30        ensembl            0 ENSG00000012048              17   human
#> 31        ensembl            0 ENSG00000012048              17   human
#>                    biotype strand                logic_name      end
#> 1           protein_coding     -1 ensembl_havana_transcript 41277387
#> 2           protein_coding     -1 ensembl_havana_transcript 41277468
#> 3           protein_coding     -1    havana_homo_sapiens_37 41277346
#> 4           protein_coding     -1    havana_homo_sapiens_37 41277346
#> 5           protein_coding     -1    havana_homo_sapiens_37 41277346
#> 6           protein_coding     -1 ensembl_havana_transcript 41277419
#> 7           protein_coding     -1 ensembl_havana_transcript 41277500
#> 8  nonsense_mediated_decay     -1    havana_homo_sapiens_37 41277305
#> 9           protein_coding     -1    havana_homo_sapiens_37 41277373
#> 10          protein_coding     -1    havana_homo_sapiens_37 41256933
#> 11          protein_coding     -1    havana_homo_sapiens_37 41277376
#> 12          protein_coding     -1    havana_homo_sapiens_37 41277419
#> 13         retained_intron     -1    havana_homo_sapiens_37 41223083
#> 14          protein_coding     -1    havana_homo_sapiens_37 41256933
#> 15          protein_coding     -1    havana_homo_sapiens_37 41243841
#> 16         retained_intron     -1    havana_homo_sapiens_37 41277332
#> 17          non_stop_decay     -1    havana_homo_sapiens_37 41247883
#> 18          protein_coding     -1    havana_homo_sapiens_37 41277340
#> 19          protein_coding     -1    havana_homo_sapiens_37 41277381
#> 20 nonsense_mediated_decay     -1    havana_homo_sapiens_37 41277317
#> 21          protein_coding     -1    havana_homo_sapiens_37 41277317
#> 22          protein_coding     -1    havana_homo_sapiens_37 41277467
#> 23          protein_coding     -1    havana_homo_sapiens_37 41256908
#> 24          protein_coding     -1    havana_homo_sapiens_37 41277370
#> 25 nonsense_mediated_decay     -1    havana_homo_sapiens_37 41277387
#> 26          protein_coding     -1    havana_homo_sapiens_37 41277338
#> 27          protein_coding     -1   ensembl_homo_sapiens_37 41277500
#> 28          protein_coding     -1   ensembl_homo_sapiens_37 41277500
#> 29          protein_coding     -1   ensembl_homo_sapiens_37 41277468
#> 30          protein_coding     -1   ensembl_homo_sapiens_37 41276132
#> 31          protein_coding     -1   ensembl_homo_sapiens_37 41277468
#>            Exon object_type Translation.db_type Translation.object_type
#> 1  c("Exon"....  Transcript                core             Translation
#> 2  c(412772....  Transcript                core             Translation
#> 3  c("core"....  Transcript                core             Translation
#> 4  c(-1, -1....  Transcript                core             Translation
#> 5  c("human....  Transcript                core             Translation
#> 6  c("core"....  Transcript                core             Translation
#> 7  c("Exon"....  Transcript                core             Translation
#> 8  c("core"....  Transcript                core             Translation
#> 9  c("Exon"....  Transcript                core             Translation
#> 10 c(412568....  Transcript                core             Translation
#> 11 c(-1, -1....  Transcript                core             Translation
#> 12 c(412772....  Transcript                core             Translation
#> 13 c("Exon"....  Transcript                <NA>                    <NA>
#> 14 c(412568....  Transcript                core             Translation
#> 15 c("human....  Transcript                core             Translation
#> 16 c(412772....  Transcript                <NA>                    <NA>
#> 17 c("ENSE0....  Transcript                core             Translation
#> 18 c("GRCh3....  Transcript                core             Translation
#> 19 c("core"....  Transcript                core             Translation
#> 20 c("core"....  Transcript                core             Translation
#> 21 c("core"....  Transcript                core             Translation
#> 22 c("Exon"....  Transcript                core             Translation
#> 23 c("human....  Transcript                core             Translation
#> 24 c("core"....  Transcript                core             Translation
#> 25 c(-1, -1....  Transcript                core             Translation
#> 26 c(412771....  Transcript                core             Translation
#> 27 c("Exon"....  Transcript                core             Translation
#> 28 c("Exon"....  Transcript                core             Translation
#> 29 c(412772....  Transcript                core             Translation
#> 30 c("Exon"....  Transcript                core             Translation
#> 31 c(-1, -1....  Transcript                core             Translation
#>    Translation.end  Translation.id Translation.species Translation.version
#> 1         41276113 ENSP00000350283               human                   3
#> 2         41276113 ENSP00000417148               human                   1
#> 3         41277202 ENSP00000465818               human                   1
#> 4         41226495 ENSP00000467329               human                   1
#> 5         41202109 ENSP00000465347               human                   1
#> 6         41258543 ENSP00000418775               human                   1
#> 7         41276113 ENSP00000418960               human                   2
#> 8         41276113 ENSP00000418548               human                   1
#> 9         41276113 ENSP00000420705               human                   2
#> 10        41256933 ENSP00000419481               human                   1
#> 11        41276113 ENSP00000420412               human                   1
#> 12        41258543 ENSP00000418819               human                   1
#> 13              NA            <NA>                <NA>                  NA
#> 14        41256933 ENSP00000418212               human                   1
#> 15        41243841 ENSP00000417241               human                   1
#> 16              NA            <NA>                <NA>                  NA
#> 17        41247883 ENSP00000397145               human                   3
#> 18        41276113 ENSP00000419274               human                   1
#> 19        41276113 ENSP00000419988               human                   1
#> 20        41276113 ENSP00000420253               human                   1
#> 21        41246659 ENSP00000418986               human                   1
#> 22        41276113 ENSP00000419103               human                   1
#> 23        41256908 ENSP00000420201               human                   1
#> 24        41276113 ENSP00000417554               human                   1
#> 25        41276113 ENSP00000417988               human                   1
#> 26        41276113 ENSP00000420781               human                   1
#> 27        41276113 ENSP00000326002               human                   6
#> 28        41276113 ENSP00000312236               human                   5
#> 29        41276113 ENSP00000246907               human                   4
#> 30        41276113 ENSP00000338007               human                   3
#> 31        41246659 ENSP00000310938               human                   4
#>    Translation.Parent Translation.length Translation.start db_type    start
#> 1     ENST00000357654               1863          41197695    core 41196312
#> 2     ENST00000468300                699          41197801    core 41196822
#> 3     ENST00000586385                173          41197695    core 41197580
#> 4     ENST00000591534                354          41197695    core 41197580
#> 5     ENST00000591849                 96          41197695    core 41197580
#> 6     ENST00000493795               1816          41197695    core 41197646
#> 7     ENST00000471181               1884          41197695    core 41197646
#> 8     ENST00000461221                 63          41256972    core 41197695
#> 9     ENST00000491747                759          41197695    core 41197695
#> 10    ENST00000484087                498          41215361    core 41215361
#> 11    ENST00000478531                623          41215361    core 41215361
#> 12    ENST00000493919                572          41215377    core 41215377
#> 13               <NA>                 NA                NA    core 41219291
#> 14    ENST00000487825                266          41228505    core 41228505
#> 15    ENST00000461574                242          41228554    core 41228554
#> 16               <NA>                 NA                NA    core 41243115
#> 17    ENST00000412061                437          41245587    core 41245587
#> 18    ENST00000470026                649          41245601    core 41245601
#> 19    ENST00000477152                622          41245603    core 41245603
#> 20    ENST00000492859                 59          41262552    core 41246129
#> 21    ENST00000497488                177          41246129    core 41246129
#> 22    ENST00000494123                473          41246129    core 41246129
#> 23    ENST00000473961                319          41246187    core 41246187
#> 24    ENST00000476777                222          41247863    core 41247863
#> 25    ENST00000461798                 63          41256972    core 41251848
#> 26    ENST00000489037                 98          41256206    core 41256206
#> 27    ENST00000354071               1598          41197695    core 41196313
#> 28    ENST00000352993                721          41197695    core 41196313
#> 29    ENST00000346315               1624          41197695    core 41196313
#> 30    ENST00000351666                680          41197695    core 41196313
#> 31    ENST00000309486               1567          41197695    core 41196313
#>    length gencode_primary version              id assembly_name display_name
#> 1    7094               0       3 ENST00000357654        GRCh37    BRCA1-001
#> 2    3273               0       1 ENST00000468300        GRCh37    BRCA1-007
#> 3     781               0       1 ENST00000586385        GRCh37    BRCA1-023
#> 4    1282               0       1 ENST00000591534        GRCh37    BRCA1-024
#> 5     563               0       1 ENST00000591849        GRCh37    BRCA1-025
#> 6    5732               0       1 ENST00000493795        GRCh37    BRCA1-006
#> 7    5936               0       2 ENST00000471181        GRCh37    BRCA1-005
#> 8    5693               0       1 ENST00000461221        GRCh37    BRCA1-010
#> 9    2379               0       2 ENST00000491747        GRCh37    BRCA1-014
#> 10   1495               0       1 ENST00000484087        GRCh37    BRCA1-015
#> 11   1972               0       1 ENST00000478531        GRCh37    BRCA1-009
#> 12   1948               0       1 ENST00000493919        GRCh37    BRCA1-008
#> 13    561               0       1 ENST00000472490        GRCh37    BRCA1-021
#> 14    800               0       1 ENST00000487825        GRCh37    BRCA1-019
#> 15    726               0       1 ENST00000461574        GRCh37    BRCA1-022
#> 16   4497               0       1 ENST00000467274        GRCh37    BRCA1-012
#> 17   1312               0       3 ENST00000412061        GRCh37    BRCA1-026
#> 18   2108               0       1 ENST00000470026        GRCh37    BRCA1-011
#> 19   1980               0       1 ENST00000477152        GRCh37    BRCA1-004
#> 20   1584               0       1 ENST00000492859        GRCh37    BRCA1-002
#> 21    779               0       1 ENST00000497488        GRCh37    BRCA1-003
#> 22   1612               0       1 ENST00000494123        GRCh37    BRCA1-013
#> 23    958               0       1 ENST00000473961        GRCh37    BRCA1-018
#> 24    769               0       1 ENST00000476777        GRCh37    BRCA1-017
#> 25    582               0       1 ENST00000461798        GRCh37    BRCA1-020
#> 26    455               0       1 ENST00000489037        GRCh37    BRCA1-016
#> 27   6411               0       3 ENST00000354071        GRCh37    BRCA1-205
#> 28   3780               0       3 ENST00000352993        GRCh37    BRCA1-204
#> 29   6451               0       3 ENST00000346315        GRCh37    BRCA1-202
#> 30   3444               0       3 ENST00000351666        GRCh37    BRCA1-203
#> 31   7114               0       4 ENST00000309486        GRCh37    BRCA1-201
#> 
#> $version
#> [1] 15
#> 
#> $start
#> [1] 41196312
#> 
#> $db_type
#> [1] "core"
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
#> $strand
#> [1] -1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
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
#> $source
#> [1] "ensembl_havana"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
#> $start
#> [1] 32889611
#> 
#> $Transcript
#>           source is_canonical          Parent      species
#> 1 ensembl_havana            0 ENSG00000139618 homo_sapiens
#> 2         havana            0 ENSG00000139618 homo_sapiens
#> 3         havana            0 ENSG00000139618 homo_sapiens
#> 4         havana            0 ENSG00000139618 homo_sapiens
#> 5         havana            0 ENSG00000139618 homo_sapiens
#> 6        ensembl            1 ENSG00000139618 homo_sapiens
#>                   biotype seq_region_name                logic_name strand
#> 1          protein_coding              13 ensembl_havana_transcript      1
#> 2          protein_coding              13    havana_homo_sapiens_37      1
#> 3 nonsense_mediated_decay              13    havana_homo_sapiens_37      1
#> 4 nonsense_mediated_decay              13    havana_homo_sapiens_37      1
#> 5         retained_intron              13    havana_homo_sapiens_37      1
#> 6          protein_coding              13   ensembl_homo_sapiens_37      1
#>        end object_type         Exon  Translation.id Translation.end
#> 1 32973347  Transcript c(328896.... ENSP00000369497        32972907
#> 2 32907428  Transcript c("core".... ENSP00000435699        32907428
#> 3 32953632  Transcript c(1, 1, .... ENSP00000433168        32950807
#> 4 32972409  Transcript c("Exon".... ENSP00000434898        32970229
#> 5 32972585  Transcript c(329709....            <NA>              NA
#> 6 32973805  Transcript c("core".... ENSP00000439902        32972907
#>   Translation.Parent Translation.length Translation.start Translation.species
#> 1    ENST00000380152               3418          32890598        homo_sapiens
#> 2    ENST00000530893                481          32899266        homo_sapiens
#> 3    ENST00000528762                 64          32945108        homo_sapiens
#> 4    ENST00000470094                186          32953977        homo_sapiens
#> 5               <NA>                 NA                NA                <NA>
#> 6    ENST00000544455               3418          32890598        homo_sapiens
#>   Translation.version Translation.db_type Translation.object_type db_type
#> 1                   3                core             Translation    core
#> 2                   2                core             Translation    core
#> 3                   1                core             Translation    core
#> 4                   1                core             Translation    core
#> 5                  NA                <NA>                    <NA>    core
#> 6                   1                core             Translation    core
#>      start length version gencode_primary              id display_name
#> 1 32889611  10930       3               0 ENST00000380152    BRCA2-001
#> 2 32889642   2011       2               0 ENST00000530893    BRCA2-003
#> 3 32945108    495       1               0 ENST00000528762    BRCA2-005
#> 4 32953977    842       1               0 ENST00000470094    BRCA2-002
#> 5 32970946    523       1               0 ENST00000533776    BRCA2-006
#> 6 32889617  10984       1               0 ENST00000544455    BRCA2-201
#>   assembly_name
#> 1        GRCh37
#> 2        GRCh37
#> 3        GRCh37
#> 4        GRCh37
#> 5        GRCh37
#> 6        GRCh37
#> 
#> $version
#> [1] 10
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $assembly_name
#> [1] "GRCh37"
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
#> $seq_region_name
#> [1] "13"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $strand
#> [1] 1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $end
#> [1] 32973805
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $end
#> [1] 41277500
#> 
#> $db_type
#> [1] "core"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $Transcript
#>    gencode_primary    start object_type strand display_name is_canonical length
#> 1                0 41196312  Transcript     -1    BRCA1-001            0   7094
#> 2                0 41196822  Transcript     -1    BRCA1-007            0   3273
#> 3                0 41197580  Transcript     -1    BRCA1-023            0    781
#> 4                0 41197580  Transcript     -1    BRCA1-024            0   1282
#> 5                0 41197580  Transcript     -1    BRCA1-025            0    563
#> 6                0 41197646  Transcript     -1    BRCA1-006            0   5732
#> 7                0 41197646  Transcript     -1    BRCA1-005            1   5936
#> 8                0 41197695  Transcript     -1    BRCA1-010            0   5693
#> 9                0 41197695  Transcript     -1    BRCA1-014            0   2379
#> 10               0 41215361  Transcript     -1    BRCA1-015            0   1495
#> 11               0 41215361  Transcript     -1    BRCA1-009            0   1972
#> 12               0 41215377  Transcript     -1    BRCA1-008            0   1948
#> 13               0 41219291  Transcript     -1    BRCA1-021            0    561
#> 14               0 41228505  Transcript     -1    BRCA1-019            0    800
#> 15               0 41228554  Transcript     -1    BRCA1-022            0    726
#> 16               0 41243115  Transcript     -1    BRCA1-012            0   4497
#> 17               0 41245587  Transcript     -1    BRCA1-026            0   1312
#> 18               0 41245601  Transcript     -1    BRCA1-011            0   2108
#> 19               0 41245603  Transcript     -1    BRCA1-004            0   1980
#> 20               0 41246129  Transcript     -1    BRCA1-002            0   1584
#> 21               0 41246129  Transcript     -1    BRCA1-003            0    779
#> 22               0 41246129  Transcript     -1    BRCA1-013            0   1612
#> 23               0 41246187  Transcript     -1    BRCA1-018            0    958
#> 24               0 41247863  Transcript     -1    BRCA1-017            0    769
#> 25               0 41251848  Transcript     -1    BRCA1-020            0    582
#> 26               0 41256206  Transcript     -1    BRCA1-016            0    455
#> 27               0 41196313  Transcript     -1    BRCA1-205            0   6411
#> 28               0 41196313  Transcript     -1    BRCA1-204            0   3780
#> 29               0 41196313  Transcript     -1    BRCA1-202            0   6451
#> 30               0 41196313  Transcript     -1    BRCA1-203            0   3444
#> 31               0 41196313  Transcript     -1    BRCA1-201            0   7114
#>    Translation.version Translation.length Translation.species  Translation.id
#> 1                    3               1863               human ENSP00000350283
#> 2                    1                699               human ENSP00000417148
#> 3                    1                173               human ENSP00000465818
#> 4                    1                354               human ENSP00000467329
#> 5                    1                 96               human ENSP00000465347
#> 6                    1               1816               human ENSP00000418775
#> 7                    2               1884               human ENSP00000418960
#> 8                    1                 63               human ENSP00000418548
#> 9                    2                759               human ENSP00000420705
#> 10                   1                498               human ENSP00000419481
#> 11                   1                623               human ENSP00000420412
#> 12                   1                572               human ENSP00000418819
#> 13                  NA                 NA                <NA>            <NA>
#> 14                   1                266               human ENSP00000418212
#> 15                   1                242               human ENSP00000417241
#> 16                  NA                 NA                <NA>            <NA>
#> 17                   3                437               human ENSP00000397145
#> 18                   1                649               human ENSP00000419274
#> 19                   1                622               human ENSP00000419988
#> 20                   1                 59               human ENSP00000420253
#> 21                   1                177               human ENSP00000418986
#> 22                   1                473               human ENSP00000419103
#> 23                   1                319               human ENSP00000420201
#> 24                   1                222               human ENSP00000417554
#> 25                   1                 63               human ENSP00000417988
#> 26                   1                 98               human ENSP00000420781
#> 27                   6               1598               human ENSP00000326002
#> 28                   5                721               human ENSP00000312236
#> 29                   4               1624               human ENSP00000246907
#> 30                   3                680               human ENSP00000338007
#> 31                   4               1567               human ENSP00000310938
#>    Translation.Parent Translation.db_type Translation.start Translation.end
#> 1     ENST00000357654                core          41197695        41276113
#> 2     ENST00000468300                core          41197801        41276113
#> 3     ENST00000586385                core          41197695        41277202
#> 4     ENST00000591534                core          41197695        41226495
#> 5     ENST00000591849                core          41197695        41202109
#> 6     ENST00000493795                core          41197695        41258543
#> 7     ENST00000471181                core          41197695        41276113
#> 8     ENST00000461221                core          41256972        41276113
#> 9     ENST00000491747                core          41197695        41276113
#> 10    ENST00000484087                core          41215361        41256933
#> 11    ENST00000478531                core          41215361        41276113
#> 12    ENST00000493919                core          41215377        41258543
#> 13               <NA>                <NA>                NA              NA
#> 14    ENST00000487825                core          41228505        41256933
#> 15    ENST00000461574                core          41228554        41243841
#> 16               <NA>                <NA>                NA              NA
#> 17    ENST00000412061                core          41245587        41247883
#> 18    ENST00000470026                core          41245601        41276113
#> 19    ENST00000477152                core          41245603        41276113
#> 20    ENST00000492859                core          41262552        41276113
#> 21    ENST00000497488                core          41246129        41246659
#> 22    ENST00000494123                core          41246129        41276113
#> 23    ENST00000473961                core          41246187        41256908
#> 24    ENST00000476777                core          41247863        41276113
#> 25    ENST00000461798                core          41256972        41276113
#> 26    ENST00000489037                core          41256206        41276113
#> 27    ENST00000354071                core          41197695        41276113
#> 28    ENST00000352993                core          41197695        41276113
#> 29    ENST00000346315                core          41197695        41276113
#> 30    ENST00000351666                core          41197695        41276113
#> 31    ENST00000309486                core          41197695        41246659
#>    Translation.object_type          Parent      end db_type         source
#> 1              Translation ENSG00000012048 41277387    core ensembl_havana
#> 2              Translation ENSG00000012048 41277468    core ensembl_havana
#> 3              Translation ENSG00000012048 41277346    core         havana
#> 4              Translation ENSG00000012048 41277346    core         havana
#> 5              Translation ENSG00000012048 41277346    core         havana
#> 6              Translation ENSG00000012048 41277419    core ensembl_havana
#> 7              Translation ENSG00000012048 41277500    core ensembl_havana
#> 8              Translation ENSG00000012048 41277305    core         havana
#> 9              Translation ENSG00000012048 41277373    core         havana
#> 10             Translation ENSG00000012048 41256933    core         havana
#> 11             Translation ENSG00000012048 41277376    core         havana
#> 12             Translation ENSG00000012048 41277419    core         havana
#> 13                    <NA> ENSG00000012048 41223083    core         havana
#> 14             Translation ENSG00000012048 41256933    core         havana
#> 15             Translation ENSG00000012048 41243841    core         havana
#> 16                    <NA> ENSG00000012048 41277332    core         havana
#> 17             Translation ENSG00000012048 41247883    core         havana
#> 18             Translation ENSG00000012048 41277340    core         havana
#> 19             Translation ENSG00000012048 41277381    core         havana
#> 20             Translation ENSG00000012048 41277317    core         havana
#> 21             Translation ENSG00000012048 41277317    core         havana
#> 22             Translation ENSG00000012048 41277467    core         havana
#> 23             Translation ENSG00000012048 41256908    core         havana
#> 24             Translation ENSG00000012048 41277370    core         havana
#> 25             Translation ENSG00000012048 41277387    core         havana
#> 26             Translation ENSG00000012048 41277338    core         havana
#> 27             Translation ENSG00000012048 41277500    core        ensembl
#> 28             Translation ENSG00000012048 41277500    core        ensembl
#> 29             Translation ENSG00000012048 41277468    core        ensembl
#> 30             Translation ENSG00000012048 41276132    core        ensembl
#> 31             Translation ENSG00000012048 41277468    core        ensembl
#>                   logic_name              id species         Exon
#> 1  ensembl_havana_transcript ENST00000357654   human c("ENSE0....
#> 2  ensembl_havana_transcript ENST00000468300   human c("human....
#> 3     havana_homo_sapiens_37 ENST00000586385   human c("human....
#> 4     havana_homo_sapiens_37 ENST00000591534   human c("GRCh3....
#> 5     havana_homo_sapiens_37 ENST00000591849   human c(-1, -1....
#> 6  ensembl_havana_transcript ENST00000493795   human c("ENSE0....
#> 7  ensembl_havana_transcript ENST00000471181   human c("GRCh3....
#> 8     havana_homo_sapiens_37 ENST00000461221   human c("ENSE0....
#> 9     havana_homo_sapiens_37 ENST00000491747   human c("human....
#> 10    havana_homo_sapiens_37 ENST00000484087   human c("ENSE0....
#> 11    havana_homo_sapiens_37 ENST00000478531   human c("Exon"....
#> 12    havana_homo_sapiens_37 ENST00000493919   human c("Exon"....
#> 13    havana_homo_sapiens_37 ENST00000472490   human c("Exon"....
#> 14    havana_homo_sapiens_37 ENST00000487825   human c("17", ....
#> 15    havana_homo_sapiens_37 ENST00000461574   human c(-1, -1....
#> 16    havana_homo_sapiens_37 ENST00000467274   human c(-1, -1....
#> 17    havana_homo_sapiens_37 ENST00000412061   human c("17", ....
#> 18    havana_homo_sapiens_37 ENST00000470026   human c("Exon"....
#> 19    havana_homo_sapiens_37 ENST00000477152   human c("ENSE0....
#> 20    havana_homo_sapiens_37 ENST00000492859   human c("Exon"....
#> 21    havana_homo_sapiens_37 ENST00000497488   human c(1, 1),....
#> 22    havana_homo_sapiens_37 ENST00000494123   human c(-1, -1....
#> 23    havana_homo_sapiens_37 ENST00000473961   human c("human....
#> 24    havana_homo_sapiens_37 ENST00000476777   human c("Exon"....
#> 25    havana_homo_sapiens_37 ENST00000461798   human c("human....
#> 26    havana_homo_sapiens_37 ENST00000489037   human c("human....
#> 27   ensembl_homo_sapiens_37 ENST00000354071   human c("human....
#> 28   ensembl_homo_sapiens_37 ENST00000352993   human c("GRCh3....
#> 29   ensembl_homo_sapiens_37 ENST00000346315   human c("Exon"....
#> 30   ensembl_homo_sapiens_37 ENST00000351666   human c("human....
#> 31   ensembl_homo_sapiens_37 ENST00000309486   human c(1, 1, ....
#>    seq_region_name version                 biotype assembly_name
#> 1               17       3          protein_coding        GRCh37
#> 2               17       1          protein_coding        GRCh37
#> 3               17       1          protein_coding        GRCh37
#> 4               17       1          protein_coding        GRCh37
#> 5               17       1          protein_coding        GRCh37
#> 6               17       1          protein_coding        GRCh37
#> 7               17       2          protein_coding        GRCh37
#> 8               17       1 nonsense_mediated_decay        GRCh37
#> 9               17       2          protein_coding        GRCh37
#> 10              17       1          protein_coding        GRCh37
#> 11              17       1          protein_coding        GRCh37
#> 12              17       1          protein_coding        GRCh37
#> 13              17       1         retained_intron        GRCh37
#> 14              17       1          protein_coding        GRCh37
#> 15              17       1          protein_coding        GRCh37
#> 16              17       1         retained_intron        GRCh37
#> 17              17       3          non_stop_decay        GRCh37
#> 18              17       1          protein_coding        GRCh37
#> 19              17       1          protein_coding        GRCh37
#> 20              17       1 nonsense_mediated_decay        GRCh37
#> 21              17       1          protein_coding        GRCh37
#> 22              17       1          protein_coding        GRCh37
#> 23              17       1          protein_coding        GRCh37
#> 24              17       1          protein_coding        GRCh37
#> 25              17       1 nonsense_mediated_decay        GRCh37
#> 26              17       1          protein_coding        GRCh37
#> 27              17       3          protein_coding        GRCh37
#> 28              17       3          protein_coding        GRCh37
#> 29              17       3          protein_coding        GRCh37
#> 30              17       3          protein_coding        GRCh37
#> 31              17       4          protein_coding        GRCh37
#> 
#> $species
#> [1] "human"
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $version
#> [1] 15
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $start
#> [1] 41196312
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $strand
#> [1] -1
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
```
