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
#> $strand
#> [1] -1
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $species
#> [1] "human"
#> 
#> $Transcript
#>    gencode_primary              id                 biotype
#> 1                0 ENST00000357654          protein_coding
#> 2                0 ENST00000468300          protein_coding
#> 3                0 ENST00000586385          protein_coding
#> 4                0 ENST00000591534          protein_coding
#> 5                0 ENST00000591849          protein_coding
#> 6                0 ENST00000493795          protein_coding
#> 7                0 ENST00000471181          protein_coding
#> 8                0 ENST00000461221 nonsense_mediated_decay
#> 9                0 ENST00000491747          protein_coding
#> 10               0 ENST00000484087          protein_coding
#> 11               0 ENST00000478531          protein_coding
#> 12               0 ENST00000493919          protein_coding
#> 13               0 ENST00000472490         retained_intron
#> 14               0 ENST00000487825          protein_coding
#> 15               0 ENST00000461574          protein_coding
#> 16               0 ENST00000467274         retained_intron
#> 17               0 ENST00000412061          non_stop_decay
#> 18               0 ENST00000470026          protein_coding
#> 19               0 ENST00000477152          protein_coding
#> 20               0 ENST00000492859 nonsense_mediated_decay
#> 21               0 ENST00000497488          protein_coding
#> 22               0 ENST00000494123          protein_coding
#> 23               0 ENST00000473961          protein_coding
#> 24               0 ENST00000476777          protein_coding
#> 25               0 ENST00000461798 nonsense_mediated_decay
#> 26               0 ENST00000489037          protein_coding
#> 27               0 ENST00000354071          protein_coding
#> 28               0 ENST00000352993          protein_coding
#> 29               0 ENST00000346315          protein_coding
#> 30               0 ENST00000351666          protein_coding
#> 31               0 ENST00000309486          protein_coding
#>                   logic_name strand Translation.Parent  Translation.id
#> 1  ensembl_havana_transcript     -1    ENST00000357654 ENSP00000350283
#> 2  ensembl_havana_transcript     -1    ENST00000468300 ENSP00000417148
#> 3     havana_homo_sapiens_37     -1    ENST00000586385 ENSP00000465818
#> 4     havana_homo_sapiens_37     -1    ENST00000591534 ENSP00000467329
#> 5     havana_homo_sapiens_37     -1    ENST00000591849 ENSP00000465347
#> 6  ensembl_havana_transcript     -1    ENST00000493795 ENSP00000418775
#> 7  ensembl_havana_transcript     -1    ENST00000471181 ENSP00000418960
#> 8     havana_homo_sapiens_37     -1    ENST00000461221 ENSP00000418548
#> 9     havana_homo_sapiens_37     -1    ENST00000491747 ENSP00000420705
#> 10    havana_homo_sapiens_37     -1    ENST00000484087 ENSP00000419481
#> 11    havana_homo_sapiens_37     -1    ENST00000478531 ENSP00000420412
#> 12    havana_homo_sapiens_37     -1    ENST00000493919 ENSP00000418819
#> 13    havana_homo_sapiens_37     -1               <NA>            <NA>
#> 14    havana_homo_sapiens_37     -1    ENST00000487825 ENSP00000418212
#> 15    havana_homo_sapiens_37     -1    ENST00000461574 ENSP00000417241
#> 16    havana_homo_sapiens_37     -1               <NA>            <NA>
#> 17    havana_homo_sapiens_37     -1    ENST00000412061 ENSP00000397145
#> 18    havana_homo_sapiens_37     -1    ENST00000470026 ENSP00000419274
#> 19    havana_homo_sapiens_37     -1    ENST00000477152 ENSP00000419988
#> 20    havana_homo_sapiens_37     -1    ENST00000492859 ENSP00000420253
#> 21    havana_homo_sapiens_37     -1    ENST00000497488 ENSP00000418986
#> 22    havana_homo_sapiens_37     -1    ENST00000494123 ENSP00000419103
#> 23    havana_homo_sapiens_37     -1    ENST00000473961 ENSP00000420201
#> 24    havana_homo_sapiens_37     -1    ENST00000476777 ENSP00000417554
#> 25    havana_homo_sapiens_37     -1    ENST00000461798 ENSP00000417988
#> 26    havana_homo_sapiens_37     -1    ENST00000489037 ENSP00000420781
#> 27   ensembl_homo_sapiens_37     -1    ENST00000354071 ENSP00000326002
#> 28   ensembl_homo_sapiens_37     -1    ENST00000352993 ENSP00000312236
#> 29   ensembl_homo_sapiens_37     -1    ENST00000346315 ENSP00000246907
#> 30   ensembl_homo_sapiens_37     -1    ENST00000351666 ENSP00000338007
#> 31   ensembl_homo_sapiens_37     -1    ENST00000309486 ENSP00000310938
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
#>    Translation.start Translation.length Translation.end Translation.version
#> 1           41197695               1863        41276113                   3
#> 2           41197801                699        41276113                   1
#> 3           41197695                173        41277202                   1
#> 4           41197695                354        41226495                   1
#> 5           41197695                 96        41202109                   1
#> 6           41197695               1816        41258543                   1
#> 7           41197695               1884        41276113                   2
#> 8           41256972                 63        41276113                   1
#> 9           41197695                759        41276113                   2
#> 10          41215361                498        41256933                   1
#> 11          41215361                623        41276113                   1
#> 12          41215377                572        41258543                   1
#> 13                NA                 NA              NA                  NA
#> 14          41228505                266        41256933                   1
#> 15          41228554                242        41243841                   1
#> 16                NA                 NA              NA                  NA
#> 17          41245587                437        41247883                   3
#> 18          41245601                649        41276113                   1
#> 19          41245603                622        41276113                   1
#> 20          41262552                 59        41276113                   1
#> 21          41246129                177        41246659                   1
#> 22          41246129                473        41276113                   1
#> 23          41246187                319        41256908                   1
#> 24          41247863                222        41276113                   1
#> 25          41256972                 63        41276113                   1
#> 26          41256206                 98        41276113                   1
#> 27          41197695               1598        41276113                   6
#> 28          41197695                721        41276113                   5
#> 29          41197695               1624        41276113                   4
#> 30          41197695                680        41276113                   3
#> 31          41197695               1567        41246659                   4
#>    is_canonical seq_region_name         Exon species         source object_type
#> 1             0              17 c("ENSE0....   human ensembl_havana  Transcript
#> 2             0              17 c(1, 1, ....   human ensembl_havana  Transcript
#> 3             0              17 c("ENSE0....   human         havana  Transcript
#> 4             0              17 c(412772....   human         havana  Transcript
#> 5             0              17 c(1, 1, ....   human         havana  Transcript
#> 6             0              17 c(1, 1, ....   human ensembl_havana  Transcript
#> 7             1              17 c("ENSE0....   human ensembl_havana  Transcript
#> 8             0              17 c(1, 1, ....   human         havana  Transcript
#> 9             0              17 c("ENSE0....   human         havana  Transcript
#> 10            0              17 c("core"....   human         havana  Transcript
#> 11            0              17 c(412772....   human         havana  Transcript
#> 12            0              17 c(1, 1, ....   human         havana  Transcript
#> 13            0              17 c(-1, -1....   human         havana  Transcript
#> 14            0              17 c("17", ....   human         havana  Transcript
#> 15            0              17 c("core"....   human         havana  Transcript
#> 16            0              17 c(1, 1, ....   human         havana  Transcript
#> 17            0              17 c("ENSE0....   human         havana  Transcript
#> 18            0              17 c("17", ....   human         havana  Transcript
#> 19            0              17 c("ENSE0....   human         havana  Transcript
#> 20            0              17 c("core"....   human         havana  Transcript
#> 21            0              17 c("ENSE0....   human         havana  Transcript
#> 22            0              17 c("ENSE0....   human         havana  Transcript
#> 23            0              17 c("ENSE0....   human         havana  Transcript
#> 24            0              17 c(412772....   human         havana  Transcript
#> 25            0              17 c("17", ....   human         havana  Transcript
#> 26            0              17 c("Exon"....   human         havana  Transcript
#> 27            0              17 c("ENSE0....   human        ensembl  Transcript
#> 28            0              17 c(412772....   human        ensembl  Transcript
#> 29            0              17 c("17", ....   human        ensembl  Transcript
#> 30            0              17 c("core"....   human        ensembl  Transcript
#> 31            0              17 c("17", ....   human        ensembl  Transcript
#>    db_type          Parent display_name length      end version    start
#> 1     core ENSG00000012048    BRCA1-001   7094 41277387       3 41196312
#> 2     core ENSG00000012048    BRCA1-007   3273 41277468       1 41196822
#> 3     core ENSG00000012048    BRCA1-023    781 41277346       1 41197580
#> 4     core ENSG00000012048    BRCA1-024   1282 41277346       1 41197580
#> 5     core ENSG00000012048    BRCA1-025    563 41277346       1 41197580
#> 6     core ENSG00000012048    BRCA1-006   5732 41277419       1 41197646
#> 7     core ENSG00000012048    BRCA1-005   5936 41277500       2 41197646
#> 8     core ENSG00000012048    BRCA1-010   5693 41277305       1 41197695
#> 9     core ENSG00000012048    BRCA1-014   2379 41277373       2 41197695
#> 10    core ENSG00000012048    BRCA1-015   1495 41256933       1 41215361
#> 11    core ENSG00000012048    BRCA1-009   1972 41277376       1 41215361
#> 12    core ENSG00000012048    BRCA1-008   1948 41277419       1 41215377
#> 13    core ENSG00000012048    BRCA1-021    561 41223083       1 41219291
#> 14    core ENSG00000012048    BRCA1-019    800 41256933       1 41228505
#> 15    core ENSG00000012048    BRCA1-022    726 41243841       1 41228554
#> 16    core ENSG00000012048    BRCA1-012   4497 41277332       1 41243115
#> 17    core ENSG00000012048    BRCA1-026   1312 41247883       3 41245587
#> 18    core ENSG00000012048    BRCA1-011   2108 41277340       1 41245601
#> 19    core ENSG00000012048    BRCA1-004   1980 41277381       1 41245603
#> 20    core ENSG00000012048    BRCA1-002   1584 41277317       1 41246129
#> 21    core ENSG00000012048    BRCA1-003    779 41277317       1 41246129
#> 22    core ENSG00000012048    BRCA1-013   1612 41277467       1 41246129
#> 23    core ENSG00000012048    BRCA1-018    958 41256908       1 41246187
#> 24    core ENSG00000012048    BRCA1-017    769 41277370       1 41247863
#> 25    core ENSG00000012048    BRCA1-020    582 41277387       1 41251848
#> 26    core ENSG00000012048    BRCA1-016    455 41277338       1 41256206
#> 27    core ENSG00000012048    BRCA1-205   6411 41277500       3 41196313
#> 28    core ENSG00000012048    BRCA1-204   3780 41277500       3 41196313
#> 29    core ENSG00000012048    BRCA1-202   6451 41277468       3 41196313
#> 30    core ENSG00000012048    BRCA1-203   3444 41276132       3 41196313
#> 31    core ENSG00000012048    BRCA1-201   7114 41277468       4 41196313
#>    assembly_name
#> 1         GRCh37
#> 2         GRCh37
#> 3         GRCh37
#> 4         GRCh37
#> 5         GRCh37
#> 6         GRCh37
#> 7         GRCh37
#> 8         GRCh37
#> 9         GRCh37
#> 10        GRCh37
#> 11        GRCh37
#> 12        GRCh37
#> 13        GRCh37
#> 14        GRCh37
#> 15        GRCh37
#> 16        GRCh37
#> 17        GRCh37
#> 18        GRCh37
#> 19        GRCh37
#> 20        GRCh37
#> 21        GRCh37
#> 22        GRCh37
#> 23        GRCh37
#> 24        GRCh37
#> 25        GRCh37
#> 26        GRCh37
#> 27        GRCh37
#> 28        GRCh37
#> 29        GRCh37
#> 30        GRCh37
#> 31        GRCh37
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $end
#> [1] 41277500
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $version
#> [1] 15
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $start
#> [1] 41196312
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $db_type
#> [1] "core"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $biotype
#> [1] "protein_coding"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $Transcript
#>            Parent object_type db_type         source assembly_name    start
#> 1 ENSG00000139618  Transcript    core ensembl_havana        GRCh37 32889611
#> 2 ENSG00000139618  Transcript    core         havana        GRCh37 32889642
#> 3 ENSG00000139618  Transcript    core         havana        GRCh37 32945108
#> 4 ENSG00000139618  Transcript    core         havana        GRCh37 32953977
#> 5 ENSG00000139618  Transcript    core         havana        GRCh37 32970946
#> 6 ENSG00000139618  Transcript    core        ensembl        GRCh37 32889617
#>   version display_name      end length              id
#> 1       3    BRCA2-001 32973347  10930 ENST00000380152
#> 2       2    BRCA2-003 32907428   2011 ENST00000530893
#> 3       1    BRCA2-005 32953632    495 ENST00000528762
#> 4       1    BRCA2-002 32972409    842 ENST00000470094
#> 5       1    BRCA2-006 32972585    523 ENST00000533776
#> 6       1    BRCA2-201 32973805  10984 ENST00000544455
#>                  logic_name                 biotype gencode_primary
#> 1 ensembl_havana_transcript          protein_coding               0
#> 2    havana_homo_sapiens_37          protein_coding               0
#> 3    havana_homo_sapiens_37 nonsense_mediated_decay               0
#> 4    havana_homo_sapiens_37 nonsense_mediated_decay               0
#> 5    havana_homo_sapiens_37         retained_intron               0
#> 6   ensembl_homo_sapiens_37          protein_coding               0
#>        species seq_region_name is_canonical         Exon strand
#> 1 homo_sapiens              13            0 c("ENSE0....      1
#> 2 homo_sapiens              13            0 c("ENSE0....      1
#> 3 homo_sapiens              13            0 c("core"....      1
#> 4 homo_sapiens              13            0 c(1, 1, ....      1
#> 5 homo_sapiens              13            0 c(329709....      1
#> 6 homo_sapiens              13            1 c("Exon"....      1
#>   Translation.version Translation.length Translation.end Translation.start
#> 1                   3               3418        32972907          32890598
#> 2                   2                481        32907428          32899266
#> 3                   1                 64        32950807          32945108
#> 4                   1                186        32970229          32953977
#> 5                  NA                 NA              NA                NA
#> 6                   1               3418        32972907          32890598
#>   Translation.species Translation.object_type Translation.db_type
#> 1        homo_sapiens             Translation                core
#> 2        homo_sapiens             Translation                core
#> 3        homo_sapiens             Translation                core
#> 4        homo_sapiens             Translation                core
#> 5                <NA>                    <NA>                <NA>
#> 6        homo_sapiens             Translation                core
#>   Translation.Parent  Translation.id
#> 1    ENST00000380152 ENSP00000369497
#> 2    ENST00000530893 ENSP00000435699
#> 3    ENST00000528762 ENSP00000433168
#> 4    ENST00000470094 ENSP00000434898
#> 5               <NA>            <NA>
#> 6    ENST00000544455 ENSP00000439902
#> 
#> $strand
#> [1] 1
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $db_type
#> [1] "core"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $start
#> [1] 32889611
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $end
#> [1] 32973805
#> 
#> $display_name
#> [1] "BRCA2"
#> 
#> $version
#> [1] 10
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
#> $start
#> [1] 41196312
#> 
#> $version
#> [1] 15
#> 
#> $Transcript
#>    is_canonical         source                logic_name strand      end
#> 1             0 ensembl_havana ensembl_havana_transcript     -1 41277387
#> 2             0 ensembl_havana ensembl_havana_transcript     -1 41277468
#> 3             0         havana    havana_homo_sapiens_37     -1 41277346
#> 4             0         havana    havana_homo_sapiens_37     -1 41277346
#> 5             0         havana    havana_homo_sapiens_37     -1 41277346
#> 6             0 ensembl_havana ensembl_havana_transcript     -1 41277419
#> 7             1 ensembl_havana ensembl_havana_transcript     -1 41277500
#> 8             0         havana    havana_homo_sapiens_37     -1 41277305
#> 9             0         havana    havana_homo_sapiens_37     -1 41277373
#> 10            0         havana    havana_homo_sapiens_37     -1 41256933
#> 11            0         havana    havana_homo_sapiens_37     -1 41277376
#> 12            0         havana    havana_homo_sapiens_37     -1 41277419
#> 13            0         havana    havana_homo_sapiens_37     -1 41223083
#> 14            0         havana    havana_homo_sapiens_37     -1 41256933
#> 15            0         havana    havana_homo_sapiens_37     -1 41243841
#> 16            0         havana    havana_homo_sapiens_37     -1 41277332
#> 17            0         havana    havana_homo_sapiens_37     -1 41247883
#> 18            0         havana    havana_homo_sapiens_37     -1 41277340
#> 19            0         havana    havana_homo_sapiens_37     -1 41277381
#> 20            0         havana    havana_homo_sapiens_37     -1 41277317
#> 21            0         havana    havana_homo_sapiens_37     -1 41277317
#> 22            0         havana    havana_homo_sapiens_37     -1 41277467
#> 23            0         havana    havana_homo_sapiens_37     -1 41256908
#> 24            0         havana    havana_homo_sapiens_37     -1 41277370
#> 25            0         havana    havana_homo_sapiens_37     -1 41277387
#> 26            0         havana    havana_homo_sapiens_37     -1 41277338
#> 27            0        ensembl   ensembl_homo_sapiens_37     -1 41277500
#> 28            0        ensembl   ensembl_homo_sapiens_37     -1 41277500
#> 29            0        ensembl   ensembl_homo_sapiens_37     -1 41277468
#> 30            0        ensembl   ensembl_homo_sapiens_37     -1 41276132
#> 31            0        ensembl   ensembl_homo_sapiens_37     -1 41277468
#>             Parent                 biotype species seq_region_name db_type
#> 1  ENSG00000012048          protein_coding   human              17    core
#> 2  ENSG00000012048          protein_coding   human              17    core
#> 3  ENSG00000012048          protein_coding   human              17    core
#> 4  ENSG00000012048          protein_coding   human              17    core
#> 5  ENSG00000012048          protein_coding   human              17    core
#> 6  ENSG00000012048          protein_coding   human              17    core
#> 7  ENSG00000012048          protein_coding   human              17    core
#> 8  ENSG00000012048 nonsense_mediated_decay   human              17    core
#> 9  ENSG00000012048          protein_coding   human              17    core
#> 10 ENSG00000012048          protein_coding   human              17    core
#> 11 ENSG00000012048          protein_coding   human              17    core
#> 12 ENSG00000012048          protein_coding   human              17    core
#> 13 ENSG00000012048         retained_intron   human              17    core
#> 14 ENSG00000012048          protein_coding   human              17    core
#> 15 ENSG00000012048          protein_coding   human              17    core
#> 16 ENSG00000012048         retained_intron   human              17    core
#> 17 ENSG00000012048          non_stop_decay   human              17    core
#> 18 ENSG00000012048          protein_coding   human              17    core
#> 19 ENSG00000012048          protein_coding   human              17    core
#> 20 ENSG00000012048 nonsense_mediated_decay   human              17    core
#> 21 ENSG00000012048          protein_coding   human              17    core
#> 22 ENSG00000012048          protein_coding   human              17    core
#> 23 ENSG00000012048          protein_coding   human              17    core
#> 24 ENSG00000012048          protein_coding   human              17    core
#> 25 ENSG00000012048 nonsense_mediated_decay   human              17    core
#> 26 ENSG00000012048          protein_coding   human              17    core
#> 27 ENSG00000012048          protein_coding   human              17    core
#> 28 ENSG00000012048          protein_coding   human              17    core
#> 29 ENSG00000012048          protein_coding   human              17    core
#> 30 ENSG00000012048          protein_coding   human              17    core
#> 31 ENSG00000012048          protein_coding   human              17    core
#>    object_type         Exon Translation.object_type Translation.db_type
#> 1   Transcript c("core"....             Translation                core
#> 2   Transcript c("core"....             Translation                core
#> 3   Transcript c("Exon"....             Translation                core
#> 4   Transcript c("17", ....             Translation                core
#> 5   Transcript c("Exon"....             Translation                core
#> 6   Transcript c("Exon"....             Translation                core
#> 7   Transcript c(1, 1, ....             Translation                core
#> 8   Transcript c("core"....             Translation                core
#> 9   Transcript c(1, 1, ....             Translation                core
#> 10  Transcript c("core"....             Translation                core
#> 11  Transcript c(412773....             Translation                core
#> 12  Transcript c("Exon"....             Translation                core
#> 13  Transcript c("Exon"....                    <NA>                <NA>
#> 14  Transcript c("Exon"....             Translation                core
#> 15  Transcript c(1, 1, ....             Translation                core
#> 16  Transcript c("core"....                    <NA>                <NA>
#> 17  Transcript c(412478....             Translation                core
#> 18  Transcript c("Exon"....             Translation                core
#> 19  Transcript c(412772....             Translation                core
#> 20  Transcript c("Exon"....             Translation                core
#> 21  Transcript c(-1, -1....             Translation                core
#> 22  Transcript c("Exon"....             Translation                core
#> 23  Transcript c(-1, -1....             Translation                core
#> 24  Transcript c("core"....             Translation                core
#> 25  Transcript c("Exon"....             Translation                core
#> 26  Transcript c(412773....             Translation                core
#> 27  Transcript c("GRCh3....             Translation                core
#> 28  Transcript c(412772....             Translation                core
#> 29  Transcript c(-1, -1....             Translation                core
#> 30  Transcript c("core"....             Translation                core
#> 31  Transcript c("ENSE0....             Translation                core
#>    Translation.species Translation.version Translation.Parent
#> 1                human                   3    ENST00000357654
#> 2                human                   1    ENST00000468300
#> 3                human                   1    ENST00000586385
#> 4                human                   1    ENST00000591534
#> 5                human                   1    ENST00000591849
#> 6                human                   1    ENST00000493795
#> 7                human                   2    ENST00000471181
#> 8                human                   1    ENST00000461221
#> 9                human                   2    ENST00000491747
#> 10               human                   1    ENST00000484087
#> 11               human                   1    ENST00000478531
#> 12               human                   1    ENST00000493919
#> 13                <NA>                  NA               <NA>
#> 14               human                   1    ENST00000487825
#> 15               human                   1    ENST00000461574
#> 16                <NA>                  NA               <NA>
#> 17               human                   3    ENST00000412061
#> 18               human                   1    ENST00000470026
#> 19               human                   1    ENST00000477152
#> 20               human                   1    ENST00000492859
#> 21               human                   1    ENST00000497488
#> 22               human                   1    ENST00000494123
#> 23               human                   1    ENST00000473961
#> 24               human                   1    ENST00000476777
#> 25               human                   1    ENST00000461798
#> 26               human                   1    ENST00000489037
#> 27               human                   6    ENST00000354071
#> 28               human                   5    ENST00000352993
#> 29               human                   4    ENST00000346315
#> 30               human                   3    ENST00000351666
#> 31               human                   4    ENST00000309486
#>    Translation.length Translation.start Translation.end  Translation.id
#> 1                1863          41197695        41276113 ENSP00000350283
#> 2                 699          41197801        41276113 ENSP00000417148
#> 3                 173          41197695        41277202 ENSP00000465818
#> 4                 354          41197695        41226495 ENSP00000467329
#> 5                  96          41197695        41202109 ENSP00000465347
#> 6                1816          41197695        41258543 ENSP00000418775
#> 7                1884          41197695        41276113 ENSP00000418960
#> 8                  63          41256972        41276113 ENSP00000418548
#> 9                 759          41197695        41276113 ENSP00000420705
#> 10                498          41215361        41256933 ENSP00000419481
#> 11                623          41215361        41276113 ENSP00000420412
#> 12                572          41215377        41258543 ENSP00000418819
#> 13                 NA                NA              NA            <NA>
#> 14                266          41228505        41256933 ENSP00000418212
#> 15                242          41228554        41243841 ENSP00000417241
#> 16                 NA                NA              NA            <NA>
#> 17                437          41245587        41247883 ENSP00000397145
#> 18                649          41245601        41276113 ENSP00000419274
#> 19                622          41245603        41276113 ENSP00000419988
#> 20                 59          41262552        41276113 ENSP00000420253
#> 21                177          41246129        41246659 ENSP00000418986
#> 22                473          41246129        41276113 ENSP00000419103
#> 23                319          41246187        41256908 ENSP00000420201
#> 24                222          41247863        41276113 ENSP00000417554
#> 25                 63          41256972        41276113 ENSP00000417988
#> 26                 98          41256206        41276113 ENSP00000420781
#> 27               1598          41197695        41276113 ENSP00000326002
#> 28                721          41197695        41276113 ENSP00000312236
#> 29               1624          41197695        41276113 ENSP00000246907
#> 30                680          41197695        41276113 ENSP00000338007
#> 31               1567          41197695        41246659 ENSP00000310938
#>                 id display_name assembly_name length    start version
#> 1  ENST00000357654    BRCA1-001        GRCh37   7094 41196312       3
#> 2  ENST00000468300    BRCA1-007        GRCh37   3273 41196822       1
#> 3  ENST00000586385    BRCA1-023        GRCh37    781 41197580       1
#> 4  ENST00000591534    BRCA1-024        GRCh37   1282 41197580       1
#> 5  ENST00000591849    BRCA1-025        GRCh37    563 41197580       1
#> 6  ENST00000493795    BRCA1-006        GRCh37   5732 41197646       1
#> 7  ENST00000471181    BRCA1-005        GRCh37   5936 41197646       2
#> 8  ENST00000461221    BRCA1-010        GRCh37   5693 41197695       1
#> 9  ENST00000491747    BRCA1-014        GRCh37   2379 41197695       2
#> 10 ENST00000484087    BRCA1-015        GRCh37   1495 41215361       1
#> 11 ENST00000478531    BRCA1-009        GRCh37   1972 41215361       1
#> 12 ENST00000493919    BRCA1-008        GRCh37   1948 41215377       1
#> 13 ENST00000472490    BRCA1-021        GRCh37    561 41219291       1
#> 14 ENST00000487825    BRCA1-019        GRCh37    800 41228505       1
#> 15 ENST00000461574    BRCA1-022        GRCh37    726 41228554       1
#> 16 ENST00000467274    BRCA1-012        GRCh37   4497 41243115       1
#> 17 ENST00000412061    BRCA1-026        GRCh37   1312 41245587       3
#> 18 ENST00000470026    BRCA1-011        GRCh37   2108 41245601       1
#> 19 ENST00000477152    BRCA1-004        GRCh37   1980 41245603       1
#> 20 ENST00000492859    BRCA1-002        GRCh37   1584 41246129       1
#> 21 ENST00000497488    BRCA1-003        GRCh37    779 41246129       1
#> 22 ENST00000494123    BRCA1-013        GRCh37   1612 41246129       1
#> 23 ENST00000473961    BRCA1-018        GRCh37    958 41246187       1
#> 24 ENST00000476777    BRCA1-017        GRCh37    769 41247863       1
#> 25 ENST00000461798    BRCA1-020        GRCh37    582 41251848       1
#> 26 ENST00000489037    BRCA1-016        GRCh37    455 41256206       1
#> 27 ENST00000354071    BRCA1-205        GRCh37   6411 41196313       3
#> 28 ENST00000352993    BRCA1-204        GRCh37   3780 41196313       3
#> 29 ENST00000346315    BRCA1-202        GRCh37   6451 41196313       3
#> 30 ENST00000351666    BRCA1-203        GRCh37   3444 41196313       3
#> 31 ENST00000309486    BRCA1-201        GRCh37   7114 41196313       4
#>    gencode_primary
#> 1                0
#> 2                0
#> 3                0
#> 4                0
#> 5                0
#> 6                0
#> 7                0
#> 8                0
#> 9                0
#> 10               0
#> 11               0
#> 12               0
#> 13               0
#> 14               0
#> 15               0
#> 16               0
#> 17               0
#> 18               0
#> 19               0
#> 20               0
#> 21               0
#> 22               0
#> 23               0
#> 24               0
#> 25               0
#> 26               0
#> 27               0
#> 28               0
#> 29               0
#> 30               0
#> 31               0
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
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
#> $strand
#> [1] -1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $end
#> [1] 41277500
#> 
```
