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
#> $id
#> [1] "ENSG00000012048"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $strand
#> [1] -1
#> 
#> $Transcript
#>    length      end display_name version    start assembly_name         source
#> 1    7094 41277387    BRCA1-001       3 41196312        GRCh37 ensembl_havana
#> 2    3273 41277468    BRCA1-007       1 41196822        GRCh37 ensembl_havana
#> 3     781 41277346    BRCA1-023       1 41197580        GRCh37         havana
#> 4    1282 41277346    BRCA1-024       1 41197580        GRCh37         havana
#> 5     563 41277346    BRCA1-025       1 41197580        GRCh37         havana
#> 6    5732 41277419    BRCA1-006       1 41197646        GRCh37 ensembl_havana
#> 7    5936 41277500    BRCA1-005       2 41197646        GRCh37 ensembl_havana
#> 8    5693 41277305    BRCA1-010       1 41197695        GRCh37         havana
#> 9    2379 41277373    BRCA1-014       2 41197695        GRCh37         havana
#> 10   1495 41256933    BRCA1-015       1 41215361        GRCh37         havana
#> 11   1972 41277376    BRCA1-009       1 41215361        GRCh37         havana
#> 12   1948 41277419    BRCA1-008       1 41215377        GRCh37         havana
#> 13    561 41223083    BRCA1-021       1 41219291        GRCh37         havana
#> 14    800 41256933    BRCA1-019       1 41228505        GRCh37         havana
#> 15    726 41243841    BRCA1-022       1 41228554        GRCh37         havana
#> 16   4497 41277332    BRCA1-012       1 41243115        GRCh37         havana
#> 17   1312 41247883    BRCA1-026       3 41245587        GRCh37         havana
#> 18   2108 41277340    BRCA1-011       1 41245601        GRCh37         havana
#> 19   1980 41277381    BRCA1-004       1 41245603        GRCh37         havana
#> 20   1584 41277317    BRCA1-002       1 41246129        GRCh37         havana
#> 21    779 41277317    BRCA1-003       1 41246129        GRCh37         havana
#> 22   1612 41277467    BRCA1-013       1 41246129        GRCh37         havana
#> 23    958 41256908    BRCA1-018       1 41246187        GRCh37         havana
#> 24    769 41277370    BRCA1-017       1 41247863        GRCh37         havana
#> 25    582 41277387    BRCA1-020       1 41251848        GRCh37         havana
#> 26    455 41277338    BRCA1-016       1 41256206        GRCh37         havana
#> 27   6411 41277500    BRCA1-205       3 41196313        GRCh37        ensembl
#> 28   3780 41277500    BRCA1-204       3 41196313        GRCh37        ensembl
#> 29   6451 41277468    BRCA1-202       3 41196313        GRCh37        ensembl
#> 30   3444 41276132    BRCA1-203       3 41196313        GRCh37        ensembl
#> 31   7114 41277468    BRCA1-201       4 41196313        GRCh37        ensembl
#>    db_type object_type          Parent Translation.Parent  Translation.id
#> 1     core  Transcript ENSG00000012048    ENST00000357654 ENSP00000350283
#> 2     core  Transcript ENSG00000012048    ENST00000468300 ENSP00000417148
#> 3     core  Transcript ENSG00000012048    ENST00000586385 ENSP00000465818
#> 4     core  Transcript ENSG00000012048    ENST00000591534 ENSP00000467329
#> 5     core  Transcript ENSG00000012048    ENST00000591849 ENSP00000465347
#> 6     core  Transcript ENSG00000012048    ENST00000493795 ENSP00000418775
#> 7     core  Transcript ENSG00000012048    ENST00000471181 ENSP00000418960
#> 8     core  Transcript ENSG00000012048    ENST00000461221 ENSP00000418548
#> 9     core  Transcript ENSG00000012048    ENST00000491747 ENSP00000420705
#> 10    core  Transcript ENSG00000012048    ENST00000484087 ENSP00000419481
#> 11    core  Transcript ENSG00000012048    ENST00000478531 ENSP00000420412
#> 12    core  Transcript ENSG00000012048    ENST00000493919 ENSP00000418819
#> 13    core  Transcript ENSG00000012048               <NA>            <NA>
#> 14    core  Transcript ENSG00000012048    ENST00000487825 ENSP00000418212
#> 15    core  Transcript ENSG00000012048    ENST00000461574 ENSP00000417241
#> 16    core  Transcript ENSG00000012048               <NA>            <NA>
#> 17    core  Transcript ENSG00000012048    ENST00000412061 ENSP00000397145
#> 18    core  Transcript ENSG00000012048    ENST00000470026 ENSP00000419274
#> 19    core  Transcript ENSG00000012048    ENST00000477152 ENSP00000419988
#> 20    core  Transcript ENSG00000012048    ENST00000492859 ENSP00000420253
#> 21    core  Transcript ENSG00000012048    ENST00000497488 ENSP00000418986
#> 22    core  Transcript ENSG00000012048    ENST00000494123 ENSP00000419103
#> 23    core  Transcript ENSG00000012048    ENST00000473961 ENSP00000420201
#> 24    core  Transcript ENSG00000012048    ENST00000476777 ENSP00000417554
#> 25    core  Transcript ENSG00000012048    ENST00000461798 ENSP00000417988
#> 26    core  Transcript ENSG00000012048    ENST00000489037 ENSP00000420781
#> 27    core  Transcript ENSG00000012048    ENST00000354071 ENSP00000326002
#> 28    core  Transcript ENSG00000012048    ENST00000352993 ENSP00000312236
#> 29    core  Transcript ENSG00000012048    ENST00000346315 ENSP00000246907
#> 30    core  Transcript ENSG00000012048    ENST00000351666 ENSP00000338007
#> 31    core  Transcript ENSG00000012048    ENST00000309486 ENSP00000310938
#>    Translation.db_type Translation.object_type Translation.species
#> 1                 core             Translation               human
#> 2                 core             Translation               human
#> 3                 core             Translation               human
#> 4                 core             Translation               human
#> 5                 core             Translation               human
#> 6                 core             Translation               human
#> 7                 core             Translation               human
#> 8                 core             Translation               human
#> 9                 core             Translation               human
#> 10                core             Translation               human
#> 11                core             Translation               human
#> 12                core             Translation               human
#> 13                <NA>                    <NA>                <NA>
#> 14                core             Translation               human
#> 15                core             Translation               human
#> 16                <NA>                    <NA>                <NA>
#> 17                core             Translation               human
#> 18                core             Translation               human
#> 19                core             Translation               human
#> 20                core             Translation               human
#> 21                core             Translation               human
#> 22                core             Translation               human
#> 23                core             Translation               human
#> 24                core             Translation               human
#> 25                core             Translation               human
#> 26                core             Translation               human
#> 27                core             Translation               human
#> 28                core             Translation               human
#> 29                core             Translation               human
#> 30                core             Translation               human
#> 31                core             Translation               human
#>    Translation.start Translation.version Translation.length Translation.end
#> 1           41197695                   3               1863        41276113
#> 2           41197801                   1                699        41276113
#> 3           41197695                   1                173        41277202
#> 4           41197695                   1                354        41226495
#> 5           41197695                   1                 96        41202109
#> 6           41197695                   1               1816        41258543
#> 7           41197695                   2               1884        41276113
#> 8           41256972                   1                 63        41276113
#> 9           41197695                   2                759        41276113
#> 10          41215361                   1                498        41256933
#> 11          41215361                   1                623        41276113
#> 12          41215377                   1                572        41258543
#> 13                NA                  NA                 NA              NA
#> 14          41228505                   1                266        41256933
#> 15          41228554                   1                242        41243841
#> 16                NA                  NA                 NA              NA
#> 17          41245587                   3                437        41247883
#> 18          41245601                   1                649        41276113
#> 19          41245603                   1                622        41276113
#> 20          41262552                   1                 59        41276113
#> 21          41246129                   1                177        41246659
#> 22          41246129                   1                473        41276113
#> 23          41246187                   1                319        41256908
#> 24          41247863                   1                222        41276113
#> 25          41256972                   1                 63        41276113
#> 26          41256206                   1                 98        41276113
#> 27          41197695                   6               1598        41276113
#> 28          41197695                   5                721        41276113
#> 29          41197695                   4               1624        41276113
#> 30          41197695                   3                680        41276113
#> 31          41197695                   4               1567        41246659
#>    strand         Exon seq_region_name is_canonical species gencode_primary
#> 1      -1 c("17", ....              17            0   human               0
#> 2      -1 c(1, 1, ....              17            0   human               0
#> 3      -1 c("ENSE0....              17            0   human               0
#> 4      -1 c("ENSE0....              17            0   human               0
#> 5      -1 c("ENSE0....              17            0   human               0
#> 6      -1 c("17", ....              17            0   human               0
#> 7      -1 c("core"....              17            1   human               0
#> 8      -1 c(-1, -1....              17            0   human               0
#> 9      -1 c(-1, -1....              17            0   human               0
#> 10     -1 c("ENSE0....              17            0   human               0
#> 11     -1 c(1, 1, ....              17            0   human               0
#> 12     -1 c("core"....              17            0   human               0
#> 13     -1 c("ENSE0....              17            0   human               0
#> 14     -1 c("core"....              17            0   human               0
#> 15     -1 c("Exon"....              17            0   human               0
#> 16     -1 c(412772....              17            0   human               0
#> 17     -1 c("17", ....              17            0   human               0
#> 18     -1 c("GRCh3....              17            0   human               0
#> 19     -1 c(412773....              17            0   human               0
#> 20     -1 c("Exon"....              17            0   human               0
#> 21     -1 c(1, 1),....              17            0   human               0
#> 22     -1 c("Exon"....              17            0   human               0
#> 23     -1 c("ENSE0....              17            0   human               0
#> 24     -1 c(1, 1, ....              17            0   human               0
#> 25     -1 c("core"....              17            0   human               0
#> 26     -1 c("GRCh3....              17            0   human               0
#> 27     -1 c("core"....              17            0   human               0
#> 28     -1 c(412772....              17            0   human               0
#> 29     -1 c("GRCh3....              17            0   human               0
#> 30     -1 c("17", ....              17            0   human               0
#> 31     -1 c(412774....              17            0   human               0
#>                   logic_name                 biotype              id
#> 1  ensembl_havana_transcript          protein_coding ENST00000357654
#> 2  ensembl_havana_transcript          protein_coding ENST00000468300
#> 3     havana_homo_sapiens_37          protein_coding ENST00000586385
#> 4     havana_homo_sapiens_37          protein_coding ENST00000591534
#> 5     havana_homo_sapiens_37          protein_coding ENST00000591849
#> 6  ensembl_havana_transcript          protein_coding ENST00000493795
#> 7  ensembl_havana_transcript          protein_coding ENST00000471181
#> 8     havana_homo_sapiens_37 nonsense_mediated_decay ENST00000461221
#> 9     havana_homo_sapiens_37          protein_coding ENST00000491747
#> 10    havana_homo_sapiens_37          protein_coding ENST00000484087
#> 11    havana_homo_sapiens_37          protein_coding ENST00000478531
#> 12    havana_homo_sapiens_37          protein_coding ENST00000493919
#> 13    havana_homo_sapiens_37         retained_intron ENST00000472490
#> 14    havana_homo_sapiens_37          protein_coding ENST00000487825
#> 15    havana_homo_sapiens_37          protein_coding ENST00000461574
#> 16    havana_homo_sapiens_37         retained_intron ENST00000467274
#> 17    havana_homo_sapiens_37          non_stop_decay ENST00000412061
#> 18    havana_homo_sapiens_37          protein_coding ENST00000470026
#> 19    havana_homo_sapiens_37          protein_coding ENST00000477152
#> 20    havana_homo_sapiens_37 nonsense_mediated_decay ENST00000492859
#> 21    havana_homo_sapiens_37          protein_coding ENST00000497488
#> 22    havana_homo_sapiens_37          protein_coding ENST00000494123
#> 23    havana_homo_sapiens_37          protein_coding ENST00000473961
#> 24    havana_homo_sapiens_37          protein_coding ENST00000476777
#> 25    havana_homo_sapiens_37 nonsense_mediated_decay ENST00000461798
#> 26    havana_homo_sapiens_37          protein_coding ENST00000489037
#> 27   ensembl_homo_sapiens_37          protein_coding ENST00000354071
#> 28   ensembl_homo_sapiens_37          protein_coding ENST00000352993
#> 29   ensembl_homo_sapiens_37          protein_coding ENST00000346315
#> 30   ensembl_homo_sapiens_37          protein_coding ENST00000351666
#> 31   ensembl_homo_sapiens_37          protein_coding ENST00000309486
#> 
#> $species
#> [1] "human"
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $version
#> [1] 15
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $end
#> [1] 41277500
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $start
#> [1] 41196312
#> 
queryEnsemblByGene("ENSG00000139618")
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
#> $Transcript
#>                id assembly_name display_name    start length gencode_primary
#> 1 ENST00000380152        GRCh37    BRCA2-001 32889611  10930               0
#> 2 ENST00000530893        GRCh37    BRCA2-003 32889642   2011               0
#> 3 ENST00000528762        GRCh37    BRCA2-005 32945108    495               0
#> 4 ENST00000470094        GRCh37    BRCA2-002 32953977    842               0
#> 5 ENST00000533776        GRCh37    BRCA2-006 32970946    523               0
#> 6 ENST00000544455        GRCh37    BRCA2-201 32889617  10984               0
#>   version db_type         Exon object_type Translation.object_type
#> 1       3    core c("GRCh3....  Transcript             Translation
#> 2       2    core c("core"....  Transcript             Translation
#> 3       1    core c("Exon"....  Transcript             Translation
#> 4       1    core c("homo_....  Transcript             Translation
#> 5       1    core c("ENSE0....  Transcript                    <NA>
#> 6       1    core c("core"....  Transcript             Translation
#>   Translation.db_type Translation.species Translation.version
#> 1                core        homo_sapiens                   3
#> 2                core        homo_sapiens                   2
#> 3                core        homo_sapiens                   1
#> 4                core        homo_sapiens                   1
#> 5                <NA>                <NA>                  NA
#> 6                core        homo_sapiens                   1
#>   Translation.Parent Translation.length Translation.start Translation.end
#> 1    ENST00000380152               3418          32890598        32972907
#> 2    ENST00000530893                481          32899266        32907428
#> 3    ENST00000528762                 64          32945108        32950807
#> 4    ENST00000470094                186          32953977        32970229
#> 5               <NA>                 NA                NA              NA
#> 6    ENST00000544455               3418          32890598        32972907
#>    Translation.id strand                logic_name      end          Parent
#> 1 ENSP00000369497      1 ensembl_havana_transcript 32973347 ENSG00000139618
#> 2 ENSP00000435699      1    havana_homo_sapiens_37 32907428 ENSG00000139618
#> 3 ENSP00000433168      1    havana_homo_sapiens_37 32953632 ENSG00000139618
#> 4 ENSP00000434898      1    havana_homo_sapiens_37 32972409 ENSG00000139618
#> 5            <NA>      1    havana_homo_sapiens_37 32972585 ENSG00000139618
#> 6 ENSP00000439902      1   ensembl_homo_sapiens_37 32973805 ENSG00000139618
#>   seq_region_name                 biotype      species is_canonical
#> 1              13          protein_coding homo_sapiens            0
#> 2              13          protein_coding homo_sapiens            0
#> 3              13 nonsense_mediated_decay homo_sapiens            0
#> 4              13 nonsense_mediated_decay homo_sapiens            0
#> 5              13         retained_intron homo_sapiens            0
#> 6              13          protein_coding homo_sapiens            1
#>           source
#> 1 ensembl_havana
#> 2         havana
#> 3         havana
#> 4         havana
#> 5         havana
#> 6        ensembl
#> 
#> $version
#> [1] 10
#> 
#> $start
#> [1] 32889611
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $source
#> [1] "ensembl_havana"
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
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $end
#> [1] 32973805
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $strand
#> [1] 1
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $object_type
#> [1] "Gene"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $db_type
#> [1] "core"
#> 
#> $end
#> [1] 41277500
#> 
#> $Transcript
#>    version species strand is_canonical seq_region_name Translation.length
#> 1        3   human     -1            0              17               1863
#> 2        1   human     -1            0              17                699
#> 3        1   human     -1            0              17                173
#> 4        1   human     -1            0              17                354
#> 5        1   human     -1            0              17                 96
#> 6        1   human     -1            0              17               1816
#> 7        2   human     -1            1              17               1884
#> 8        1   human     -1            0              17                 63
#> 9        2   human     -1            0              17                759
#> 10       1   human     -1            0              17                498
#> 11       1   human     -1            0              17                623
#> 12       1   human     -1            0              17                572
#> 13       1   human     -1            0              17                 NA
#> 14       1   human     -1            0              17                266
#> 15       1   human     -1            0              17                242
#> 16       1   human     -1            0              17                 NA
#> 17       3   human     -1            0              17                437
#> 18       1   human     -1            0              17                649
#> 19       1   human     -1            0              17                622
#> 20       1   human     -1            0              17                 59
#> 21       1   human     -1            0              17                177
#> 22       1   human     -1            0              17                473
#> 23       1   human     -1            0              17                319
#> 24       1   human     -1            0              17                222
#> 25       1   human     -1            0              17                 63
#> 26       1   human     -1            0              17                 98
#> 27       3   human     -1            0              17               1598
#> 28       3   human     -1            0              17                721
#> 29       3   human     -1            0              17               1624
#> 30       3   human     -1            0              17                680
#> 31       4   human     -1            0              17               1567
#>     Translation.id Translation.Parent Translation.object_type Translation.start
#> 1  ENSP00000350283    ENST00000357654             Translation          41197695
#> 2  ENSP00000417148    ENST00000468300             Translation          41197801
#> 3  ENSP00000465818    ENST00000586385             Translation          41197695
#> 4  ENSP00000467329    ENST00000591534             Translation          41197695
#> 5  ENSP00000465347    ENST00000591849             Translation          41197695
#> 6  ENSP00000418775    ENST00000493795             Translation          41197695
#> 7  ENSP00000418960    ENST00000471181             Translation          41197695
#> 8  ENSP00000418548    ENST00000461221             Translation          41256972
#> 9  ENSP00000420705    ENST00000491747             Translation          41197695
#> 10 ENSP00000419481    ENST00000484087             Translation          41215361
#> 11 ENSP00000420412    ENST00000478531             Translation          41215361
#> 12 ENSP00000418819    ENST00000493919             Translation          41215377
#> 13            <NA>               <NA>                    <NA>                NA
#> 14 ENSP00000418212    ENST00000487825             Translation          41228505
#> 15 ENSP00000417241    ENST00000461574             Translation          41228554
#> 16            <NA>               <NA>                    <NA>                NA
#> 17 ENSP00000397145    ENST00000412061             Translation          41245587
#> 18 ENSP00000419274    ENST00000470026             Translation          41245601
#> 19 ENSP00000419988    ENST00000477152             Translation          41245603
#> 20 ENSP00000420253    ENST00000492859             Translation          41262552
#> 21 ENSP00000418986    ENST00000497488             Translation          41246129
#> 22 ENSP00000419103    ENST00000494123             Translation          41246129
#> 23 ENSP00000420201    ENST00000473961             Translation          41246187
#> 24 ENSP00000417554    ENST00000476777             Translation          41247863
#> 25 ENSP00000417988    ENST00000461798             Translation          41256972
#> 26 ENSP00000420781    ENST00000489037             Translation          41256206
#> 27 ENSP00000326002    ENST00000354071             Translation          41197695
#> 28 ENSP00000312236    ENST00000352993             Translation          41197695
#> 29 ENSP00000246907    ENST00000346315             Translation          41197695
#> 30 ENSP00000338007    ENST00000351666             Translation          41197695
#> 31 ENSP00000310938    ENST00000309486             Translation          41197695
#>    Translation.end Translation.db_type Translation.version Translation.species
#> 1         41276113                core                   3               human
#> 2         41276113                core                   1               human
#> 3         41277202                core                   1               human
#> 4         41226495                core                   1               human
#> 5         41202109                core                   1               human
#> 6         41258543                core                   1               human
#> 7         41276113                core                   2               human
#> 8         41276113                core                   1               human
#> 9         41276113                core                   2               human
#> 10        41256933                core                   1               human
#> 11        41276113                core                   1               human
#> 12        41258543                core                   1               human
#> 13              NA                <NA>                  NA                <NA>
#> 14        41256933                core                   1               human
#> 15        41243841                core                   1               human
#> 16              NA                <NA>                  NA                <NA>
#> 17        41247883                core                   3               human
#> 18        41276113                core                   1               human
#> 19        41276113                core                   1               human
#> 20        41276113                core                   1               human
#> 21        41246659                core                   1               human
#> 22        41276113                core                   1               human
#> 23        41256908                core                   1               human
#> 24        41276113                core                   1               human
#> 25        41276113                core                   1               human
#> 26        41276113                core                   1               human
#> 27        41276113                core                   6               human
#> 28        41276113                core                   5               human
#> 29        41276113                core                   4               human
#> 30        41276113                core                   3               human
#> 31        41246659                core                   4               human
#>                 id         source          Parent assembly_name display_name
#> 1  ENST00000357654 ensembl_havana ENSG00000012048        GRCh37    BRCA1-001
#> 2  ENST00000468300 ensembl_havana ENSG00000012048        GRCh37    BRCA1-007
#> 3  ENST00000586385         havana ENSG00000012048        GRCh37    BRCA1-023
#> 4  ENST00000591534         havana ENSG00000012048        GRCh37    BRCA1-024
#> 5  ENST00000591849         havana ENSG00000012048        GRCh37    BRCA1-025
#> 6  ENST00000493795 ensembl_havana ENSG00000012048        GRCh37    BRCA1-006
#> 7  ENST00000471181 ensembl_havana ENSG00000012048        GRCh37    BRCA1-005
#> 8  ENST00000461221         havana ENSG00000012048        GRCh37    BRCA1-010
#> 9  ENST00000491747         havana ENSG00000012048        GRCh37    BRCA1-014
#> 10 ENST00000484087         havana ENSG00000012048        GRCh37    BRCA1-015
#> 11 ENST00000478531         havana ENSG00000012048        GRCh37    BRCA1-009
#> 12 ENST00000493919         havana ENSG00000012048        GRCh37    BRCA1-008
#> 13 ENST00000472490         havana ENSG00000012048        GRCh37    BRCA1-021
#> 14 ENST00000487825         havana ENSG00000012048        GRCh37    BRCA1-019
#> 15 ENST00000461574         havana ENSG00000012048        GRCh37    BRCA1-022
#> 16 ENST00000467274         havana ENSG00000012048        GRCh37    BRCA1-012
#> 17 ENST00000412061         havana ENSG00000012048        GRCh37    BRCA1-026
#> 18 ENST00000470026         havana ENSG00000012048        GRCh37    BRCA1-011
#> 19 ENST00000477152         havana ENSG00000012048        GRCh37    BRCA1-004
#> 20 ENST00000492859         havana ENSG00000012048        GRCh37    BRCA1-002
#> 21 ENST00000497488         havana ENSG00000012048        GRCh37    BRCA1-003
#> 22 ENST00000494123         havana ENSG00000012048        GRCh37    BRCA1-013
#> 23 ENST00000473961         havana ENSG00000012048        GRCh37    BRCA1-018
#> 24 ENST00000476777         havana ENSG00000012048        GRCh37    BRCA1-017
#> 25 ENST00000461798         havana ENSG00000012048        GRCh37    BRCA1-020
#> 26 ENST00000489037         havana ENSG00000012048        GRCh37    BRCA1-016
#> 27 ENST00000354071        ensembl ENSG00000012048        GRCh37    BRCA1-205
#> 28 ENST00000352993        ensembl ENSG00000012048        GRCh37    BRCA1-204
#> 29 ENST00000346315        ensembl ENSG00000012048        GRCh37    BRCA1-202
#> 30 ENST00000351666        ensembl ENSG00000012048        GRCh37    BRCA1-203
#> 31 ENST00000309486        ensembl ENSG00000012048        GRCh37    BRCA1-201
#>            Exon gencode_primary                logic_name length
#> 1  c(-1, -1....               0 ensembl_havana_transcript   7094
#> 2  c(1, 1, ....               0 ensembl_havana_transcript   3273
#> 3  c(1, 1, ....               0    havana_homo_sapiens_37    781
#> 4  c(-1, -1....               0    havana_homo_sapiens_37   1282
#> 5  c("Exon"....               0    havana_homo_sapiens_37    563
#> 6  c("GRCh3....               0 ensembl_havana_transcript   5732
#> 7  c("human....               0 ensembl_havana_transcript   5936
#> 8  c("Exon"....               0    havana_homo_sapiens_37   5693
#> 9  c("GRCh3....               0    havana_homo_sapiens_37   2379
#> 10 c(1, 1, ....               0    havana_homo_sapiens_37   1495
#> 11 c("GRCh3....               0    havana_homo_sapiens_37   1972
#> 12 c(-1, -1....               0    havana_homo_sapiens_37   1948
#> 13 c(1, 1),....               0    havana_homo_sapiens_37    561
#> 14 c(-1, -1....               0    havana_homo_sapiens_37    800
#> 15 c(412434....               0    havana_homo_sapiens_37    726
#> 16 c("core"....               0    havana_homo_sapiens_37   4497
#> 17 c("human....               0    havana_homo_sapiens_37   1312
#> 18 c(412773....               0    havana_homo_sapiens_37   2108
#> 19 c(412773....               0    havana_homo_sapiens_37   1980
#> 20 c(412772....               0    havana_homo_sapiens_37   1584
#> 21 c(412773....               0    havana_homo_sapiens_37    779
#> 22 c("human....               0    havana_homo_sapiens_37   1612
#> 23 c(1, 1, ....               0    havana_homo_sapiens_37    958
#> 24 c(1, 1, ....               0    havana_homo_sapiens_37    769
#> 25 c(1, 1, ....               0    havana_homo_sapiens_37    582
#> 26 c("Exon"....               0    havana_homo_sapiens_37    455
#> 27 c(-1, -1....               0   ensembl_homo_sapiens_37   6411
#> 28 c(1, 1, ....               0   ensembl_homo_sapiens_37   3780
#> 29 c(412772....               0   ensembl_homo_sapiens_37   6451
#> 30 c("Exon"....               0   ensembl_homo_sapiens_37   3444
#> 31 c(412774....               0   ensembl_homo_sapiens_37   7114
#>                    biotype object_type    start db_type      end
#> 1           protein_coding  Transcript 41196312    core 41277387
#> 2           protein_coding  Transcript 41196822    core 41277468
#> 3           protein_coding  Transcript 41197580    core 41277346
#> 4           protein_coding  Transcript 41197580    core 41277346
#> 5           protein_coding  Transcript 41197580    core 41277346
#> 6           protein_coding  Transcript 41197646    core 41277419
#> 7           protein_coding  Transcript 41197646    core 41277500
#> 8  nonsense_mediated_decay  Transcript 41197695    core 41277305
#> 9           protein_coding  Transcript 41197695    core 41277373
#> 10          protein_coding  Transcript 41215361    core 41256933
#> 11          protein_coding  Transcript 41215361    core 41277376
#> 12          protein_coding  Transcript 41215377    core 41277419
#> 13         retained_intron  Transcript 41219291    core 41223083
#> 14          protein_coding  Transcript 41228505    core 41256933
#> 15          protein_coding  Transcript 41228554    core 41243841
#> 16         retained_intron  Transcript 41243115    core 41277332
#> 17          non_stop_decay  Transcript 41245587    core 41247883
#> 18          protein_coding  Transcript 41245601    core 41277340
#> 19          protein_coding  Transcript 41245603    core 41277381
#> 20 nonsense_mediated_decay  Transcript 41246129    core 41277317
#> 21          protein_coding  Transcript 41246129    core 41277317
#> 22          protein_coding  Transcript 41246129    core 41277467
#> 23          protein_coding  Transcript 41246187    core 41256908
#> 24          protein_coding  Transcript 41247863    core 41277370
#> 25 nonsense_mediated_decay  Transcript 41251848    core 41277387
#> 26          protein_coding  Transcript 41256206    core 41277338
#> 27          protein_coding  Transcript 41196313    core 41277500
#> 28          protein_coding  Transcript 41196313    core 41277500
#> 29          protein_coding  Transcript 41196313    core 41277468
#> 30          protein_coding  Transcript 41196313    core 41276132
#> 31          protein_coding  Transcript 41196313    core 41277468
#> 
#> $start
#> [1] 41196312
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
#> $seq_region_name
#> [1] "17"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $species
#> [1] "human"
#> 
#> $strand
#> [1] -1
#> 
#> $version
#> [1] 15
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
```
