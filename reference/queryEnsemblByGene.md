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
#> $Transcript
#>             Parent gencode_primary                 biotype Translation.end
#> 1  ENSG00000012048               0          protein_coding        41276113
#> 2  ENSG00000012048               0          protein_coding        41276113
#> 3  ENSG00000012048               0          protein_coding        41277202
#> 4  ENSG00000012048               0          protein_coding        41226495
#> 5  ENSG00000012048               0          protein_coding        41202109
#> 6  ENSG00000012048               0          protein_coding        41258543
#> 7  ENSG00000012048               0          protein_coding        41276113
#> 8  ENSG00000012048               0 nonsense_mediated_decay        41276113
#> 9  ENSG00000012048               0          protein_coding        41276113
#> 10 ENSG00000012048               0          protein_coding        41256933
#> 11 ENSG00000012048               0          protein_coding        41276113
#> 12 ENSG00000012048               0          protein_coding        41258543
#> 13 ENSG00000012048               0         retained_intron              NA
#> 14 ENSG00000012048               0          protein_coding        41256933
#> 15 ENSG00000012048               0          protein_coding        41243841
#> 16 ENSG00000012048               0         retained_intron              NA
#> 17 ENSG00000012048               0          non_stop_decay        41247883
#> 18 ENSG00000012048               0          protein_coding        41276113
#> 19 ENSG00000012048               0          protein_coding        41276113
#> 20 ENSG00000012048               0 nonsense_mediated_decay        41276113
#> 21 ENSG00000012048               0          protein_coding        41246659
#> 22 ENSG00000012048               0          protein_coding        41276113
#> 23 ENSG00000012048               0          protein_coding        41256908
#> 24 ENSG00000012048               0          protein_coding        41276113
#> 25 ENSG00000012048               0 nonsense_mediated_decay        41276113
#> 26 ENSG00000012048               0          protein_coding        41276113
#> 27 ENSG00000012048               0          protein_coding        41276113
#> 28 ENSG00000012048               0          protein_coding        41276113
#> 29 ENSG00000012048               0          protein_coding        41276113
#> 30 ENSG00000012048               0          protein_coding        41276113
#> 31 ENSG00000012048               0          protein_coding        41246659
#>    Translation.object_type Translation.species Translation.start
#> 1              Translation               human          41197695
#> 2              Translation               human          41197801
#> 3              Translation               human          41197695
#> 4              Translation               human          41197695
#> 5              Translation               human          41197695
#> 6              Translation               human          41197695
#> 7              Translation               human          41197695
#> 8              Translation               human          41256972
#> 9              Translation               human          41197695
#> 10             Translation               human          41215361
#> 11             Translation               human          41215361
#> 12             Translation               human          41215377
#> 13                    <NA>                <NA>                NA
#> 14             Translation               human          41228505
#> 15             Translation               human          41228554
#> 16                    <NA>                <NA>                NA
#> 17             Translation               human          41245587
#> 18             Translation               human          41245601
#> 19             Translation               human          41245603
#> 20             Translation               human          41262552
#> 21             Translation               human          41246129
#> 22             Translation               human          41246129
#> 23             Translation               human          41246187
#> 24             Translation               human          41247863
#> 25             Translation               human          41256972
#> 26             Translation               human          41256206
#> 27             Translation               human          41197695
#> 28             Translation               human          41197695
#> 29             Translation               human          41197695
#> 30             Translation               human          41197695
#> 31             Translation               human          41197695
#>    Translation.length Translation.db_type Translation.Parent  Translation.id
#> 1                1863                core    ENST00000357654 ENSP00000350283
#> 2                 699                core    ENST00000468300 ENSP00000417148
#> 3                 173                core    ENST00000586385 ENSP00000465818
#> 4                 354                core    ENST00000591534 ENSP00000467329
#> 5                  96                core    ENST00000591849 ENSP00000465347
#> 6                1816                core    ENST00000493795 ENSP00000418775
#> 7                1884                core    ENST00000471181 ENSP00000418960
#> 8                  63                core    ENST00000461221 ENSP00000418548
#> 9                 759                core    ENST00000491747 ENSP00000420705
#> 10                498                core    ENST00000484087 ENSP00000419481
#> 11                623                core    ENST00000478531 ENSP00000420412
#> 12                572                core    ENST00000493919 ENSP00000418819
#> 13                 NA                <NA>               <NA>            <NA>
#> 14                266                core    ENST00000487825 ENSP00000418212
#> 15                242                core    ENST00000461574 ENSP00000417241
#> 16                 NA                <NA>               <NA>            <NA>
#> 17                437                core    ENST00000412061 ENSP00000397145
#> 18                649                core    ENST00000470026 ENSP00000419274
#> 19                622                core    ENST00000477152 ENSP00000419988
#> 20                 59                core    ENST00000492859 ENSP00000420253
#> 21                177                core    ENST00000497488 ENSP00000418986
#> 22                473                core    ENST00000494123 ENSP00000419103
#> 23                319                core    ENST00000473961 ENSP00000420201
#> 24                222                core    ENST00000476777 ENSP00000417554
#> 25                 63                core    ENST00000461798 ENSP00000417988
#> 26                 98                core    ENST00000489037 ENSP00000420781
#> 27               1598                core    ENST00000354071 ENSP00000326002
#> 28                721                core    ENST00000352993 ENSP00000312236
#> 29               1624                core    ENST00000346315 ENSP00000246907
#> 30                680                core    ENST00000351666 ENSP00000338007
#> 31               1567                core    ENST00000309486 ENSP00000310938
#>    Translation.version         source length    start is_canonical
#> 1                    3 ensembl_havana   7094 41196312            0
#> 2                    1 ensembl_havana   3273 41196822            0
#> 3                    1         havana    781 41197580            0
#> 4                    1         havana   1282 41197580            0
#> 5                    1         havana    563 41197580            0
#> 6                    1 ensembl_havana   5732 41197646            0
#> 7                    2 ensembl_havana   5936 41197646            1
#> 8                    1         havana   5693 41197695            0
#> 9                    2         havana   2379 41197695            0
#> 10                   1         havana   1495 41215361            0
#> 11                   1         havana   1972 41215361            0
#> 12                   1         havana   1948 41215377            0
#> 13                  NA         havana    561 41219291            0
#> 14                   1         havana    800 41228505            0
#> 15                   1         havana    726 41228554            0
#> 16                  NA         havana   4497 41243115            0
#> 17                   3         havana   1312 41245587            0
#> 18                   1         havana   2108 41245601            0
#> 19                   1         havana   1980 41245603            0
#> 20                   1         havana   1584 41246129            0
#> 21                   1         havana    779 41246129            0
#> 22                   1         havana   1612 41246129            0
#> 23                   1         havana    958 41246187            0
#> 24                   1         havana    769 41247863            0
#> 25                   1         havana    582 41251848            0
#> 26                   1         havana    455 41256206            0
#> 27                   6        ensembl   6411 41196313            0
#> 28                   5        ensembl   3780 41196313            0
#> 29                   4        ensembl   6451 41196313            0
#> 30                   3        ensembl   3444 41196313            0
#> 31                   4        ensembl   7114 41196313            0
#>                   logic_name db_type seq_region_name strand version
#> 1  ensembl_havana_transcript    core              17     -1       3
#> 2  ensembl_havana_transcript    core              17     -1       1
#> 3     havana_homo_sapiens_37    core              17     -1       1
#> 4     havana_homo_sapiens_37    core              17     -1       1
#> 5     havana_homo_sapiens_37    core              17     -1       1
#> 6  ensembl_havana_transcript    core              17     -1       1
#> 7  ensembl_havana_transcript    core              17     -1       2
#> 8     havana_homo_sapiens_37    core              17     -1       1
#> 9     havana_homo_sapiens_37    core              17     -1       2
#> 10    havana_homo_sapiens_37    core              17     -1       1
#> 11    havana_homo_sapiens_37    core              17     -1       1
#> 12    havana_homo_sapiens_37    core              17     -1       1
#> 13    havana_homo_sapiens_37    core              17     -1       1
#> 14    havana_homo_sapiens_37    core              17     -1       1
#> 15    havana_homo_sapiens_37    core              17     -1       1
#> 16    havana_homo_sapiens_37    core              17     -1       1
#> 17    havana_homo_sapiens_37    core              17     -1       3
#> 18    havana_homo_sapiens_37    core              17     -1       1
#> 19    havana_homo_sapiens_37    core              17     -1       1
#> 20    havana_homo_sapiens_37    core              17     -1       1
#> 21    havana_homo_sapiens_37    core              17     -1       1
#> 22    havana_homo_sapiens_37    core              17     -1       1
#> 23    havana_homo_sapiens_37    core              17     -1       1
#> 24    havana_homo_sapiens_37    core              17     -1       1
#> 25    havana_homo_sapiens_37    core              17     -1       1
#> 26    havana_homo_sapiens_37    core              17     -1       1
#> 27   ensembl_homo_sapiens_37    core              17     -1       3
#> 28   ensembl_homo_sapiens_37    core              17     -1       3
#> 29   ensembl_homo_sapiens_37    core              17     -1       3
#> 30   ensembl_homo_sapiens_37    core              17     -1       3
#> 31   ensembl_homo_sapiens_37    core              17     -1       4
#>                 id assembly_name         Exon display_name      end object_type
#> 1  ENST00000357654        GRCh37 c("GRCh3....    BRCA1-001 41277387  Transcript
#> 2  ENST00000468300        GRCh37 c(1, 1, ....    BRCA1-007 41277468  Transcript
#> 3  ENST00000586385        GRCh37 c("core"....    BRCA1-023 41277346  Transcript
#> 4  ENST00000591534        GRCh37 c("human....    BRCA1-024 41277346  Transcript
#> 5  ENST00000591849        GRCh37 c("Exon"....    BRCA1-025 41277346  Transcript
#> 6  ENST00000493795        GRCh37 c("core"....    BRCA1-006 41277419  Transcript
#> 7  ENST00000471181        GRCh37 c("Exon"....    BRCA1-005 41277500  Transcript
#> 8  ENST00000461221        GRCh37 c(412773....    BRCA1-010 41277305  Transcript
#> 9  ENST00000491747        GRCh37 c("GRCh3....    BRCA1-014 41277373  Transcript
#> 10 ENST00000484087        GRCh37 c("core"....    BRCA1-015 41256933  Transcript
#> 11 ENST00000478531        GRCh37 c(412773....    BRCA1-009 41277376  Transcript
#> 12 ENST00000493919        GRCh37 c("human....    BRCA1-008 41277419  Transcript
#> 13 ENST00000472490        GRCh37 c("GRCh3....    BRCA1-021 41223083  Transcript
#> 14 ENST00000487825        GRCh37 c("GRCh3....    BRCA1-019 41256933  Transcript
#> 15 ENST00000461574        GRCh37 c(412434....    BRCA1-022 41243841  Transcript
#> 16 ENST00000467274        GRCh37 c(-1, -1....    BRCA1-012 41277332  Transcript
#> 17 ENST00000412061        GRCh37 c(-1, -1....    BRCA1-026 41247883  Transcript
#> 18 ENST00000470026        GRCh37 c("core"....    BRCA1-011 41277340  Transcript
#> 19 ENST00000477152        GRCh37 c(412773....    BRCA1-004 41277381  Transcript
#> 20 ENST00000492859        GRCh37 c(1, 1, ....    BRCA1-002 41277317  Transcript
#> 21 ENST00000497488        GRCh37 c("core"....    BRCA1-003 41277317  Transcript
#> 22 ENST00000494123        GRCh37 c("core"....    BRCA1-013 41277467  Transcript
#> 23 ENST00000473961        GRCh37 c("GRCh3....    BRCA1-018 41256908  Transcript
#> 24 ENST00000476777        GRCh37 c(412773....    BRCA1-017 41277370  Transcript
#> 25 ENST00000461798        GRCh37 c("human....    BRCA1-020 41277387  Transcript
#> 26 ENST00000489037        GRCh37 c(-1, -1....    BRCA1-016 41277338  Transcript
#> 27 ENST00000354071        GRCh37 c("GRCh3....    BRCA1-205 41277500  Transcript
#> 28 ENST00000352993        GRCh37 c(412775....    BRCA1-204 41277500  Transcript
#> 29 ENST00000346315        GRCh37 c(412774....    BRCA1-202 41277468  Transcript
#> 30 ENST00000351666        GRCh37 c(412761....    BRCA1-203 41276132  Transcript
#> 31 ENST00000309486        GRCh37 c(-1, -1....    BRCA1-201 41277468  Transcript
#>    species
#> 1    human
#> 2    human
#> 3    human
#> 4    human
#> 5    human
#> 6    human
#> 7    human
#> 8    human
#> 9    human
#> 10   human
#> 11   human
#> 12   human
#> 13   human
#> 14   human
#> 15   human
#> 16   human
#> 17   human
#> 18   human
#> 19   human
#> 20   human
#> 21   human
#> 22   human
#> 23   human
#> 24   human
#> 25   human
#> 26   human
#> 27   human
#> 28   human
#> 29   human
#> 30   human
#> 31   human
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $db_type
#> [1] "core"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $strand
#> [1] -1
#> 
#> $version
#> [1] 15
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $end
#> [1] 41277500
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $species
#> [1] "human"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $start
#> [1] 41196312
#> 
#> $source
#> [1] "ensembl_havana"
#> 
queryEnsemblByGene("ENSG00000139618")
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
#> $id
#> [1] "ENSG00000139618"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $Transcript
#>   display_name      end length version    start assembly_name         source
#> 1    BRCA2-001 32973347  10930       3 32889611        GRCh37 ensembl_havana
#> 2    BRCA2-003 32907428   2011       2 32889642        GRCh37         havana
#> 3    BRCA2-005 32953632    495       1 32945108        GRCh37         havana
#> 4    BRCA2-002 32972409    842       1 32953977        GRCh37         havana
#> 5    BRCA2-006 32972585    523       1 32970946        GRCh37         havana
#> 6    BRCA2-201 32973805  10984       1 32889617        GRCh37        ensembl
#>   object_type db_type          Parent strand Translation.length Translation.end
#> 1  Transcript    core ENSG00000139618      1               3418        32972907
#> 2  Transcript    core ENSG00000139618      1                481        32907428
#> 3  Transcript    core ENSG00000139618      1                 64        32950807
#> 4  Transcript    core ENSG00000139618      1                186        32970229
#> 5  Transcript    core ENSG00000139618      1                 NA              NA
#> 6  Transcript    core ENSG00000139618      1               3418        32972907
#>   Translation.version Translation.start Translation.species
#> 1                   3          32890598        homo_sapiens
#> 2                   2          32899266        homo_sapiens
#> 3                   1          32945108        homo_sapiens
#> 4                   1          32953977        homo_sapiens
#> 5                  NA                NA                <NA>
#> 6                   1          32890598        homo_sapiens
#>   Translation.object_type Translation.db_type Translation.Parent
#> 1             Translation                core    ENST00000380152
#> 2             Translation                core    ENST00000530893
#> 3             Translation                core    ENST00000528762
#> 4             Translation                core    ENST00000470094
#> 5                    <NA>                <NA>               <NA>
#> 6             Translation                core    ENST00000544455
#>    Translation.id is_canonical seq_region_name         Exon      species
#> 1 ENSP00000369497            0              13 c("core".... homo_sapiens
#> 2 ENSP00000435699            0              13 c(1, 1, .... homo_sapiens
#> 3 ENSP00000433168            0              13 c("ENSE0.... homo_sapiens
#> 4 ENSP00000434898            0              13 c("Exon".... homo_sapiens
#> 5            <NA>            0              13 c("ENSE0.... homo_sapiens
#> 6 ENSP00000439902            1              13 c("ENSE0.... homo_sapiens
#>   gencode_primary              id                 biotype
#> 1               0 ENST00000380152          protein_coding
#> 2               0 ENST00000530893          protein_coding
#> 3               0 ENST00000528762 nonsense_mediated_decay
#> 4               0 ENST00000470094 nonsense_mediated_decay
#> 5               0 ENST00000533776         retained_intron
#> 6               0 ENST00000544455          protein_coding
#>                  logic_name
#> 1 ensembl_havana_transcript
#> 2    havana_homo_sapiens_37
#> 3    havana_homo_sapiens_37
#> 4    havana_homo_sapiens_37
#> 5    havana_homo_sapiens_37
#> 6   ensembl_homo_sapiens_37
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $strand
#> [1] 1
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $display_name
#> [1] "BRCA1"
#> 
#> $strand
#> [1] -1
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $start
#> [1] 41196312
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $end
#> [1] 41277500
#> 
#> $db_type
#> [1] "core"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $version
#> [1] 15
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $species
#> [1] "human"
#> 
#> $Transcript
#>    display_name strand gencode_primary object_type    start length
#> 1     BRCA1-001     -1               0  Transcript 41196312   7094
#> 2     BRCA1-007     -1               0  Transcript 41196822   3273
#> 3     BRCA1-023     -1               0  Transcript 41197580    781
#> 4     BRCA1-024     -1               0  Transcript 41197580   1282
#> 5     BRCA1-025     -1               0  Transcript 41197580    563
#> 6     BRCA1-006     -1               0  Transcript 41197646   5732
#> 7     BRCA1-005     -1               0  Transcript 41197646   5936
#> 8     BRCA1-010     -1               0  Transcript 41197695   5693
#> 9     BRCA1-014     -1               0  Transcript 41197695   2379
#> 10    BRCA1-015     -1               0  Transcript 41215361   1495
#> 11    BRCA1-009     -1               0  Transcript 41215361   1972
#> 12    BRCA1-008     -1               0  Transcript 41215377   1948
#> 13    BRCA1-021     -1               0  Transcript 41219291    561
#> 14    BRCA1-019     -1               0  Transcript 41228505    800
#> 15    BRCA1-022     -1               0  Transcript 41228554    726
#> 16    BRCA1-012     -1               0  Transcript 41243115   4497
#> 17    BRCA1-026     -1               0  Transcript 41245587   1312
#> 18    BRCA1-011     -1               0  Transcript 41245601   2108
#> 19    BRCA1-004     -1               0  Transcript 41245603   1980
#> 20    BRCA1-002     -1               0  Transcript 41246129   1584
#> 21    BRCA1-003     -1               0  Transcript 41246129    779
#> 22    BRCA1-013     -1               0  Transcript 41246129   1612
#> 23    BRCA1-018     -1               0  Transcript 41246187    958
#> 24    BRCA1-017     -1               0  Transcript 41247863    769
#> 25    BRCA1-020     -1               0  Transcript 41251848    582
#> 26    BRCA1-016     -1               0  Transcript 41256206    455
#> 27    BRCA1-205     -1               0  Transcript 41196313   6411
#> 28    BRCA1-204     -1               0  Transcript 41196313   3780
#> 29    BRCA1-202     -1               0  Transcript 41196313   6451
#> 30    BRCA1-203     -1               0  Transcript 41196313   3444
#> 31    BRCA1-201     -1               0  Transcript 41196313   7114
#>    Translation.db_type Translation.Parent Translation.end Translation.start
#> 1                 core    ENST00000357654        41276113          41197695
#> 2                 core    ENST00000468300        41276113          41197801
#> 3                 core    ENST00000586385        41277202          41197695
#> 4                 core    ENST00000591534        41226495          41197695
#> 5                 core    ENST00000591849        41202109          41197695
#> 6                 core    ENST00000493795        41258543          41197695
#> 7                 core    ENST00000471181        41276113          41197695
#> 8                 core    ENST00000461221        41276113          41256972
#> 9                 core    ENST00000491747        41276113          41197695
#> 10                core    ENST00000484087        41256933          41215361
#> 11                core    ENST00000478531        41276113          41215361
#> 12                core    ENST00000493919        41258543          41215377
#> 13                <NA>               <NA>              NA                NA
#> 14                core    ENST00000487825        41256933          41228505
#> 15                core    ENST00000461574        41243841          41228554
#> 16                <NA>               <NA>              NA                NA
#> 17                core    ENST00000412061        41247883          41245587
#> 18                core    ENST00000470026        41276113          41245601
#> 19                core    ENST00000477152        41276113          41245603
#> 20                core    ENST00000492859        41276113          41262552
#> 21                core    ENST00000497488        41246659          41246129
#> 22                core    ENST00000494123        41276113          41246129
#> 23                core    ENST00000473961        41256908          41246187
#> 24                core    ENST00000476777        41276113          41247863
#> 25                core    ENST00000461798        41276113          41256972
#> 26                core    ENST00000489037        41276113          41256206
#> 27                core    ENST00000354071        41276113          41197695
#> 28                core    ENST00000352993        41276113          41197695
#> 29                core    ENST00000346315        41276113          41197695
#> 30                core    ENST00000351666        41276113          41197695
#> 31                core    ENST00000309486        41246659          41197695
#>    Translation.object_type  Translation.id Translation.species
#> 1              Translation ENSP00000350283               human
#> 2              Translation ENSP00000417148               human
#> 3              Translation ENSP00000465818               human
#> 4              Translation ENSP00000467329               human
#> 5              Translation ENSP00000465347               human
#> 6              Translation ENSP00000418775               human
#> 7              Translation ENSP00000418960               human
#> 8              Translation ENSP00000418548               human
#> 9              Translation ENSP00000420705               human
#> 10             Translation ENSP00000419481               human
#> 11             Translation ENSP00000420412               human
#> 12             Translation ENSP00000418819               human
#> 13                    <NA>            <NA>                <NA>
#> 14             Translation ENSP00000418212               human
#> 15             Translation ENSP00000417241               human
#> 16                    <NA>            <NA>                <NA>
#> 17             Translation ENSP00000397145               human
#> 18             Translation ENSP00000419274               human
#> 19             Translation ENSP00000419988               human
#> 20             Translation ENSP00000420253               human
#> 21             Translation ENSP00000418986               human
#> 22             Translation ENSP00000419103               human
#> 23             Translation ENSP00000420201               human
#> 24             Translation ENSP00000417554               human
#> 25             Translation ENSP00000417988               human
#> 26             Translation ENSP00000420781               human
#> 27             Translation ENSP00000326002               human
#> 28             Translation ENSP00000312236               human
#> 29             Translation ENSP00000246907               human
#> 30             Translation ENSP00000338007               human
#> 31             Translation ENSP00000310938               human
#>    Translation.version Translation.length is_canonical              id
#> 1                    3               1863            0 ENST00000357654
#> 2                    1                699            0 ENST00000468300
#> 3                    1                173            0 ENST00000586385
#> 4                    1                354            0 ENST00000591534
#> 5                    1                 96            0 ENST00000591849
#> 6                    1               1816            0 ENST00000493795
#> 7                    2               1884            1 ENST00000471181
#> 8                    1                 63            0 ENST00000461221
#> 9                    2                759            0 ENST00000491747
#> 10                   1                498            0 ENST00000484087
#> 11                   1                623            0 ENST00000478531
#> 12                   1                572            0 ENST00000493919
#> 13                  NA                 NA            0 ENST00000472490
#> 14                   1                266            0 ENST00000487825
#> 15                   1                242            0 ENST00000461574
#> 16                  NA                 NA            0 ENST00000467274
#> 17                   3                437            0 ENST00000412061
#> 18                   1                649            0 ENST00000470026
#> 19                   1                622            0 ENST00000477152
#> 20                   1                 59            0 ENST00000492859
#> 21                   1                177            0 ENST00000497488
#> 22                   1                473            0 ENST00000494123
#> 23                   1                319            0 ENST00000473961
#> 24                   1                222            0 ENST00000476777
#> 25                   1                 63            0 ENST00000461798
#> 26                   1                 98            0 ENST00000489037
#> 27                   6               1598            0 ENST00000354071
#> 28                   5                721            0 ENST00000352993
#> 29                   4               1624            0 ENST00000346315
#> 30                   3                680            0 ENST00000351666
#> 31                   4               1567            0 ENST00000309486
#>                   logic_name         source          Parent      end db_type
#> 1  ensembl_havana_transcript ensembl_havana ENSG00000012048 41277387    core
#> 2  ensembl_havana_transcript ensembl_havana ENSG00000012048 41277468    core
#> 3     havana_homo_sapiens_37         havana ENSG00000012048 41277346    core
#> 4     havana_homo_sapiens_37         havana ENSG00000012048 41277346    core
#> 5     havana_homo_sapiens_37         havana ENSG00000012048 41277346    core
#> 6  ensembl_havana_transcript ensembl_havana ENSG00000012048 41277419    core
#> 7  ensembl_havana_transcript ensembl_havana ENSG00000012048 41277500    core
#> 8     havana_homo_sapiens_37         havana ENSG00000012048 41277305    core
#> 9     havana_homo_sapiens_37         havana ENSG00000012048 41277373    core
#> 10    havana_homo_sapiens_37         havana ENSG00000012048 41256933    core
#> 11    havana_homo_sapiens_37         havana ENSG00000012048 41277376    core
#> 12    havana_homo_sapiens_37         havana ENSG00000012048 41277419    core
#> 13    havana_homo_sapiens_37         havana ENSG00000012048 41223083    core
#> 14    havana_homo_sapiens_37         havana ENSG00000012048 41256933    core
#> 15    havana_homo_sapiens_37         havana ENSG00000012048 41243841    core
#> 16    havana_homo_sapiens_37         havana ENSG00000012048 41277332    core
#> 17    havana_homo_sapiens_37         havana ENSG00000012048 41247883    core
#> 18    havana_homo_sapiens_37         havana ENSG00000012048 41277340    core
#> 19    havana_homo_sapiens_37         havana ENSG00000012048 41277381    core
#> 20    havana_homo_sapiens_37         havana ENSG00000012048 41277317    core
#> 21    havana_homo_sapiens_37         havana ENSG00000012048 41277317    core
#> 22    havana_homo_sapiens_37         havana ENSG00000012048 41277467    core
#> 23    havana_homo_sapiens_37         havana ENSG00000012048 41256908    core
#> 24    havana_homo_sapiens_37         havana ENSG00000012048 41277370    core
#> 25    havana_homo_sapiens_37         havana ENSG00000012048 41277387    core
#> 26    havana_homo_sapiens_37         havana ENSG00000012048 41277338    core
#> 27   ensembl_homo_sapiens_37        ensembl ENSG00000012048 41277500    core
#> 28   ensembl_homo_sapiens_37        ensembl ENSG00000012048 41277500    core
#> 29   ensembl_homo_sapiens_37        ensembl ENSG00000012048 41277468    core
#> 30   ensembl_homo_sapiens_37        ensembl ENSG00000012048 41276132    core
#> 31   ensembl_homo_sapiens_37        ensembl ENSG00000012048 41277468    core
#>    seq_region_name assembly_name                 biotype version         Exon
#> 1               17        GRCh37          protein_coding       3 c(-1, -1....
#> 2               17        GRCh37          protein_coding       1 c("human....
#> 3               17        GRCh37          protein_coding       1 c(-1, -1....
#> 4               17        GRCh37          protein_coding       1 c("Exon"....
#> 5               17        GRCh37          protein_coding       1 c("human....
#> 6               17        GRCh37          protein_coding       1 c("human....
#> 7               17        GRCh37          protein_coding       2 c("human....
#> 8               17        GRCh37 nonsense_mediated_decay       1 c(-1, -1....
#> 9               17        GRCh37          protein_coding       2 c(-1, -1....
#> 10              17        GRCh37          protein_coding       1 c(-1, -1....
#> 11              17        GRCh37          protein_coding       1 c(-1, -1....
#> 12              17        GRCh37          protein_coding       1 c(1, 1, ....
#> 13              17        GRCh37         retained_intron       1 c("human....
#> 14              17        GRCh37          protein_coding       1 c(-1, -1....
#> 15              17        GRCh37          protein_coding       1 c(-1, -1....
#> 16              17        GRCh37         retained_intron       1 c("core"....
#> 17              17        GRCh37          non_stop_decay       3 c("human....
#> 18              17        GRCh37          protein_coding       1 c("17", ....
#> 19              17        GRCh37          protein_coding       1 c("human....
#> 20              17        GRCh37 nonsense_mediated_decay       1 c(412772....
#> 21              17        GRCh37          protein_coding       1 c(412773....
#> 22              17        GRCh37          protein_coding       1 c("human....
#> 23              17        GRCh37          protein_coding       1 c(-1, -1....
#> 24              17        GRCh37          protein_coding       1 c(-1, -1....
#> 25              17        GRCh37 nonsense_mediated_decay       1 c("human....
#> 26              17        GRCh37          protein_coding       1 c("GRCh3....
#> 27              17        GRCh37          protein_coding       3 c("human....
#> 28              17        GRCh37          protein_coding       3 c("human....
#> 29              17        GRCh37          protein_coding       3 c("Exon"....
#> 30              17        GRCh37          protein_coding       3 c(-1, -1....
#> 31              17        GRCh37          protein_coding       4 c("Exon"....
#>    species
#> 1    human
#> 2    human
#> 3    human
#> 4    human
#> 5    human
#> 6    human
#> 7    human
#> 8    human
#> 9    human
#> 10   human
#> 11   human
#> 12   human
#> 13   human
#> 14   human
#> 15   human
#> 16   human
#> 17   human
#> 18   human
#> 19   human
#> 20   human
#> 21   human
#> 22   human
#> 23   human
#> 24   human
#> 25   human
#> 26   human
#> 27   human
#> 28   human
#> 29   human
#> 30   human
#> 31   human
#> 
```
