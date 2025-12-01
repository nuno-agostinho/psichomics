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
#> $seq_region_name
#> [1] "17"
#> 
#> $Transcript
#>            source db_type object_type          Parent      end length
#> 1  ensembl_havana    core  Transcript ENSG00000012048 41277387   7094
#> 2  ensembl_havana    core  Transcript ENSG00000012048 41277468   3273
#> 3          havana    core  Transcript ENSG00000012048 41277346    781
#> 4          havana    core  Transcript ENSG00000012048 41277346   1282
#> 5          havana    core  Transcript ENSG00000012048 41277346    563
#> 6  ensembl_havana    core  Transcript ENSG00000012048 41277419   5732
#> 7  ensembl_havana    core  Transcript ENSG00000012048 41277500   5936
#> 8          havana    core  Transcript ENSG00000012048 41277305   5693
#> 9          havana    core  Transcript ENSG00000012048 41277373   2379
#> 10         havana    core  Transcript ENSG00000012048 41256933   1495
#> 11         havana    core  Transcript ENSG00000012048 41277376   1972
#> 12         havana    core  Transcript ENSG00000012048 41277419   1948
#> 13         havana    core  Transcript ENSG00000012048 41223083    561
#> 14         havana    core  Transcript ENSG00000012048 41256933    800
#> 15         havana    core  Transcript ENSG00000012048 41243841    726
#> 16         havana    core  Transcript ENSG00000012048 41277332   4497
#> 17         havana    core  Transcript ENSG00000012048 41247883   1312
#> 18         havana    core  Transcript ENSG00000012048 41277340   2108
#> 19         havana    core  Transcript ENSG00000012048 41277381   1980
#> 20         havana    core  Transcript ENSG00000012048 41277317   1584
#> 21         havana    core  Transcript ENSG00000012048 41277317    779
#> 22         havana    core  Transcript ENSG00000012048 41277467   1612
#> 23         havana    core  Transcript ENSG00000012048 41256908    958
#> 24         havana    core  Transcript ENSG00000012048 41277370    769
#> 25         havana    core  Transcript ENSG00000012048 41277387    582
#> 26         havana    core  Transcript ENSG00000012048 41277338    455
#> 27        ensembl    core  Transcript ENSG00000012048 41277500   6411
#> 28        ensembl    core  Transcript ENSG00000012048 41277500   3780
#> 29        ensembl    core  Transcript ENSG00000012048 41277468   6451
#> 30        ensembl    core  Transcript ENSG00000012048 41276132   3444
#> 31        ensembl    core  Transcript ENSG00000012048 41277468   7114
#>    display_name version    start assembly_name gencode_primary
#> 1     BRCA1-001       3 41196312        GRCh37               0
#> 2     BRCA1-007       1 41196822        GRCh37               0
#> 3     BRCA1-023       1 41197580        GRCh37               0
#> 4     BRCA1-024       1 41197580        GRCh37               0
#> 5     BRCA1-025       1 41197580        GRCh37               0
#> 6     BRCA1-006       1 41197646        GRCh37               0
#> 7     BRCA1-005       2 41197646        GRCh37               0
#> 8     BRCA1-010       1 41197695        GRCh37               0
#> 9     BRCA1-014       2 41197695        GRCh37               0
#> 10    BRCA1-015       1 41215361        GRCh37               0
#> 11    BRCA1-009       1 41215361        GRCh37               0
#> 12    BRCA1-008       1 41215377        GRCh37               0
#> 13    BRCA1-021       1 41219291        GRCh37               0
#> 14    BRCA1-019       1 41228505        GRCh37               0
#> 15    BRCA1-022       1 41228554        GRCh37               0
#> 16    BRCA1-012       1 41243115        GRCh37               0
#> 17    BRCA1-026       3 41245587        GRCh37               0
#> 18    BRCA1-011       1 41245601        GRCh37               0
#> 19    BRCA1-004       1 41245603        GRCh37               0
#> 20    BRCA1-002       1 41246129        GRCh37               0
#> 21    BRCA1-003       1 41246129        GRCh37               0
#> 22    BRCA1-013       1 41246129        GRCh37               0
#> 23    BRCA1-018       1 41246187        GRCh37               0
#> 24    BRCA1-017       1 41247863        GRCh37               0
#> 25    BRCA1-020       1 41251848        GRCh37               0
#> 26    BRCA1-016       1 41256206        GRCh37               0
#> 27    BRCA1-205       3 41196313        GRCh37               0
#> 28    BRCA1-204       3 41196313        GRCh37               0
#> 29    BRCA1-202       3 41196313        GRCh37               0
#> 30    BRCA1-203       3 41196313        GRCh37               0
#> 31    BRCA1-201       4 41196313        GRCh37               0
#>                    biotype                logic_name              id
#> 1           protein_coding ensembl_havana_transcript ENST00000357654
#> 2           protein_coding ensembl_havana_transcript ENST00000468300
#> 3           protein_coding    havana_homo_sapiens_37 ENST00000586385
#> 4           protein_coding    havana_homo_sapiens_37 ENST00000591534
#> 5           protein_coding    havana_homo_sapiens_37 ENST00000591849
#> 6           protein_coding ensembl_havana_transcript ENST00000493795
#> 7           protein_coding ensembl_havana_transcript ENST00000471181
#> 8  nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000461221
#> 9           protein_coding    havana_homo_sapiens_37 ENST00000491747
#> 10          protein_coding    havana_homo_sapiens_37 ENST00000484087
#> 11          protein_coding    havana_homo_sapiens_37 ENST00000478531
#> 12          protein_coding    havana_homo_sapiens_37 ENST00000493919
#> 13         retained_intron    havana_homo_sapiens_37 ENST00000472490
#> 14          protein_coding    havana_homo_sapiens_37 ENST00000487825
#> 15          protein_coding    havana_homo_sapiens_37 ENST00000461574
#> 16         retained_intron    havana_homo_sapiens_37 ENST00000467274
#> 17          non_stop_decay    havana_homo_sapiens_37 ENST00000412061
#> 18          protein_coding    havana_homo_sapiens_37 ENST00000470026
#> 19          protein_coding    havana_homo_sapiens_37 ENST00000477152
#> 20 nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000492859
#> 21          protein_coding    havana_homo_sapiens_37 ENST00000497488
#> 22          protein_coding    havana_homo_sapiens_37 ENST00000494123
#> 23          protein_coding    havana_homo_sapiens_37 ENST00000473961
#> 24          protein_coding    havana_homo_sapiens_37 ENST00000476777
#> 25 nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000461798
#> 26          protein_coding    havana_homo_sapiens_37 ENST00000489037
#> 27          protein_coding   ensembl_homo_sapiens_37 ENST00000354071
#> 28          protein_coding   ensembl_homo_sapiens_37 ENST00000352993
#> 29          protein_coding   ensembl_homo_sapiens_37 ENST00000346315
#> 30          protein_coding   ensembl_homo_sapiens_37 ENST00000351666
#> 31          protein_coding   ensembl_homo_sapiens_37 ENST00000309486
#>    Translation.version Translation.length Translation.end Translation.species
#> 1                    3               1863        41276113               human
#> 2                    1                699        41276113               human
#> 3                    1                173        41277202               human
#> 4                    1                354        41226495               human
#> 5                    1                 96        41202109               human
#> 6                    1               1816        41258543               human
#> 7                    2               1884        41276113               human
#> 8                    1                 63        41276113               human
#> 9                    2                759        41276113               human
#> 10                   1                498        41256933               human
#> 11                   1                623        41276113               human
#> 12                   1                572        41258543               human
#> 13                  NA                 NA              NA                <NA>
#> 14                   1                266        41256933               human
#> 15                   1                242        41243841               human
#> 16                  NA                 NA              NA                <NA>
#> 17                   3                437        41247883               human
#> 18                   1                649        41276113               human
#> 19                   1                622        41276113               human
#> 20                   1                 59        41276113               human
#> 21                   1                177        41246659               human
#> 22                   1                473        41276113               human
#> 23                   1                319        41256908               human
#> 24                   1                222        41276113               human
#> 25                   1                 63        41276113               human
#> 26                   1                 98        41276113               human
#> 27                   6               1598        41276113               human
#> 28                   5                721        41276113               human
#> 29                   4               1624        41276113               human
#> 30                   3                680        41276113               human
#> 31                   4               1567        41246659               human
#>    Translation.start Translation.object_type Translation.db_type
#> 1           41197695             Translation                core
#> 2           41197801             Translation                core
#> 3           41197695             Translation                core
#> 4           41197695             Translation                core
#> 5           41197695             Translation                core
#> 6           41197695             Translation                core
#> 7           41197695             Translation                core
#> 8           41256972             Translation                core
#> 9           41197695             Translation                core
#> 10          41215361             Translation                core
#> 11          41215361             Translation                core
#> 12          41215377             Translation                core
#> 13                NA                    <NA>                <NA>
#> 14          41228505             Translation                core
#> 15          41228554             Translation                core
#> 16                NA                    <NA>                <NA>
#> 17          41245587             Translation                core
#> 18          41245601             Translation                core
#> 19          41245603             Translation                core
#> 20          41262552             Translation                core
#> 21          41246129             Translation                core
#> 22          41246129             Translation                core
#> 23          41246187             Translation                core
#> 24          41247863             Translation                core
#> 25          41256972             Translation                core
#> 26          41256206             Translation                core
#> 27          41197695             Translation                core
#> 28          41197695             Translation                core
#> 29          41197695             Translation                core
#> 30          41197695             Translation                core
#> 31          41197695             Translation                core
#>     Translation.id Translation.Parent strand         Exon is_canonical
#> 1  ENSP00000350283    ENST00000357654     -1 c("core"....            0
#> 2  ENSP00000417148    ENST00000468300     -1 c("ENSE0....            0
#> 3  ENSP00000465818    ENST00000586385     -1 c(412771....            0
#> 4  ENSP00000467329    ENST00000591534     -1 c("17", ....            0
#> 5  ENSP00000465347    ENST00000591849     -1 c("GRCh3....            0
#> 6  ENSP00000418775    ENST00000493795     -1 c("ENSE0....            0
#> 7  ENSP00000418960    ENST00000471181     -1 c("core"....            1
#> 8  ENSP00000418548    ENST00000461221     -1 c("Exon"....            0
#> 9  ENSP00000420705    ENST00000491747     -1 c("Exon"....            0
#> 10 ENSP00000419481    ENST00000484087     -1 c("human....            0
#> 11 ENSP00000420412    ENST00000478531     -1 c(412773....            0
#> 12 ENSP00000418819    ENST00000493919     -1 c("ENSE0....            0
#> 13            <NA>               <NA>     -1 c("Exon"....            0
#> 14 ENSP00000418212    ENST00000487825     -1 c("17", ....            0
#> 15 ENSP00000417241    ENST00000461574     -1 c("Exon"....            0
#> 16            <NA>               <NA>     -1 c("core"....            0
#> 17 ENSP00000397145    ENST00000412061     -1 c(1, 1),....            0
#> 18 ENSP00000419274    ENST00000470026     -1 c(1, 1, ....            0
#> 19 ENSP00000419988    ENST00000477152     -1 c(1, 1, ....            0
#> 20 ENSP00000420253    ENST00000492859     -1 c("17", ....            0
#> 21 ENSP00000418986    ENST00000497488     -1 c(412773....            0
#> 22 ENSP00000419103    ENST00000494123     -1 c("ENSE0....            0
#> 23 ENSP00000420201    ENST00000473961     -1 c("human....            0
#> 24 ENSP00000417554    ENST00000476777     -1 c(-1, -1....            0
#> 25 ENSP00000417988    ENST00000461798     -1 c("Exon"....            0
#> 26 ENSP00000420781    ENST00000489037     -1 c("ENSE0....            0
#> 27 ENSP00000326002    ENST00000354071     -1 c(1, 1, ....            0
#> 28 ENSP00000312236    ENST00000352993     -1 c("Exon"....            0
#> 29 ENSP00000246907    ENST00000346315     -1 c("ENSE0....            0
#> 30 ENSP00000338007    ENST00000351666     -1 c("17", ....            0
#> 31 ENSP00000310938    ENST00000309486     -1 c(412772....            0
#>    seq_region_name species
#> 1               17   human
#> 2               17   human
#> 3               17   human
#> 4               17   human
#> 5               17   human
#> 6               17   human
#> 7               17   human
#> 8               17   human
#> 9               17   human
#> 10              17   human
#> 11              17   human
#> 12              17   human
#> 13              17   human
#> 14              17   human
#> 15              17   human
#> 16              17   human
#> 17              17   human
#> 18              17   human
#> 19              17   human
#> 20              17   human
#> 21              17   human
#> 22              17   human
#> 23              17   human
#> 24              17   human
#> 25              17   human
#> 26              17   human
#> 27              17   human
#> 28              17   human
#> 29              17   human
#> 30              17   human
#> 31              17   human
#> 
#> $species
#> [1] "human"
#> 
#> $strand
#> [1] -1
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $start
#> [1] 41196312
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $end
#> [1] 41277500
#> 
#> $version
#> [1] 15
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
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
queryEnsemblByGene("ENSG00000139618")
#> $seq_region_name
#> [1] "13"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $species
#> [1] "homo_sapiens"
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
#> $strand
#> [1] 1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $Transcript
#>           source is_canonical                 biotype      species
#> 1 ensembl_havana            0          protein_coding homo_sapiens
#> 2         havana            0          protein_coding homo_sapiens
#> 3         havana            0 nonsense_mediated_decay homo_sapiens
#> 4         havana            0 nonsense_mediated_decay homo_sapiens
#> 5         havana            0         retained_intron homo_sapiens
#> 6        ensembl            1          protein_coding homo_sapiens
#>   seq_region_name          Parent      end                logic_name strand
#> 1              13 ENSG00000139618 32973347 ensembl_havana_transcript      1
#> 2              13 ENSG00000139618 32907428    havana_homo_sapiens_37      1
#> 3              13 ENSG00000139618 32953632    havana_homo_sapiens_37      1
#> 4              13 ENSG00000139618 32972409    havana_homo_sapiens_37      1
#> 5              13 ENSG00000139618 32972585    havana_homo_sapiens_37      1
#> 6              13 ENSG00000139618 32973805   ensembl_homo_sapiens_37      1
#>    Translation.id Translation.end Translation.length Translation.Parent
#> 1 ENSP00000369497        32972907               3418    ENST00000380152
#> 2 ENSP00000435699        32907428                481    ENST00000530893
#> 3 ENSP00000433168        32950807                 64    ENST00000528762
#> 4 ENSP00000434898        32970229                186    ENST00000470094
#> 5            <NA>              NA                 NA               <NA>
#> 6 ENSP00000439902        32972907               3418    ENST00000544455
#>   Translation.start Translation.species Translation.version Translation.db_type
#> 1          32890598        homo_sapiens                   3                core
#> 2          32899266        homo_sapiens                   2                core
#> 3          32945108        homo_sapiens                   1                core
#> 4          32953977        homo_sapiens                   1                core
#> 5                NA                <NA>                  NA                <NA>
#> 6          32890598        homo_sapiens                   1                core
#>   Translation.object_type object_type         Exon db_type version
#> 1             Translation  Transcript c("core"....    core       3
#> 2             Translation  Transcript c("Exon"....    core       2
#> 3             Translation  Transcript c("Exon"....    core       1
#> 4             Translation  Transcript c(1, 1, ....    core       1
#> 5                    <NA>  Transcript c("Exon"....    core       1
#> 6             Translation  Transcript c(328898....    core       1
#>   gencode_primary    start length display_name assembly_name              id
#> 1               0 32889611  10930    BRCA2-001        GRCh37 ENST00000380152
#> 2               0 32889642   2011    BRCA2-003        GRCh37 ENST00000530893
#> 3               0 32945108    495    BRCA2-005        GRCh37 ENST00000528762
#> 4               0 32953977    842    BRCA2-002        GRCh37 ENST00000470094
#> 5               0 32970946    523    BRCA2-006        GRCh37 ENST00000533776
#> 6               0 32889617  10984    BRCA2-201        GRCh37 ENST00000544455
#> 
#> $version
#> [1] 10
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
#> $id
#> [1] "ENSG00000139618"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $db_type
#> [1] "core"
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
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
#> $start
#> [1] 41196312
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
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
#>       start assembly_name      end length display_name version          Parent
#> 1  41196312        GRCh37 41277387   7094    BRCA1-001       3 ENSG00000012048
#> 2  41196822        GRCh37 41277468   3273    BRCA1-007       1 ENSG00000012048
#> 3  41197580        GRCh37 41277346    781    BRCA1-023       1 ENSG00000012048
#> 4  41197580        GRCh37 41277346   1282    BRCA1-024       1 ENSG00000012048
#> 5  41197580        GRCh37 41277346    563    BRCA1-025       1 ENSG00000012048
#> 6  41197646        GRCh37 41277419   5732    BRCA1-006       1 ENSG00000012048
#> 7  41197646        GRCh37 41277500   5936    BRCA1-005       2 ENSG00000012048
#> 8  41197695        GRCh37 41277305   5693    BRCA1-010       1 ENSG00000012048
#> 9  41197695        GRCh37 41277373   2379    BRCA1-014       2 ENSG00000012048
#> 10 41215361        GRCh37 41256933   1495    BRCA1-015       1 ENSG00000012048
#> 11 41215361        GRCh37 41277376   1972    BRCA1-009       1 ENSG00000012048
#> 12 41215377        GRCh37 41277419   1948    BRCA1-008       1 ENSG00000012048
#> 13 41219291        GRCh37 41223083    561    BRCA1-021       1 ENSG00000012048
#> 14 41228505        GRCh37 41256933    800    BRCA1-019       1 ENSG00000012048
#> 15 41228554        GRCh37 41243841    726    BRCA1-022       1 ENSG00000012048
#> 16 41243115        GRCh37 41277332   4497    BRCA1-012       1 ENSG00000012048
#> 17 41245587        GRCh37 41247883   1312    BRCA1-026       3 ENSG00000012048
#> 18 41245601        GRCh37 41277340   2108    BRCA1-011       1 ENSG00000012048
#> 19 41245603        GRCh37 41277381   1980    BRCA1-004       1 ENSG00000012048
#> 20 41246129        GRCh37 41277317   1584    BRCA1-002       1 ENSG00000012048
#> 21 41246129        GRCh37 41277317    779    BRCA1-003       1 ENSG00000012048
#> 22 41246129        GRCh37 41277467   1612    BRCA1-013       1 ENSG00000012048
#> 23 41246187        GRCh37 41256908    958    BRCA1-018       1 ENSG00000012048
#> 24 41247863        GRCh37 41277370    769    BRCA1-017       1 ENSG00000012048
#> 25 41251848        GRCh37 41277387    582    BRCA1-020       1 ENSG00000012048
#> 26 41256206        GRCh37 41277338    455    BRCA1-016       1 ENSG00000012048
#> 27 41196313        GRCh37 41277500   6411    BRCA1-205       3 ENSG00000012048
#> 28 41196313        GRCh37 41277500   3780    BRCA1-204       3 ENSG00000012048
#> 29 41196313        GRCh37 41277468   6451    BRCA1-202       3 ENSG00000012048
#> 30 41196313        GRCh37 41276132   3444    BRCA1-203       3 ENSG00000012048
#> 31 41196313        GRCh37 41277468   7114    BRCA1-201       4 ENSG00000012048
#>            source db_type object_type         Exon seq_region_name is_canonical
#> 1  ensembl_havana    core  Transcript c("17", ....              17            0
#> 2  ensembl_havana    core  Transcript c("ENSE0....              17            0
#> 3          havana    core  Transcript c(-1, -1....              17            0
#> 4          havana    core  Transcript c("17", ....              17            0
#> 5          havana    core  Transcript c(412772....              17            0
#> 6  ensembl_havana    core  Transcript c("core"....              17            0
#> 7  ensembl_havana    core  Transcript c(1, 1, ....              17            1
#> 8          havana    core  Transcript c("core"....              17            0
#> 9          havana    core  Transcript c(412773....              17            0
#> 10         havana    core  Transcript c("ENSE0....              17            0
#> 11         havana    core  Transcript c(412772....              17            0
#> 12         havana    core  Transcript c(412772....              17            0
#> 13         havana    core  Transcript c(1, 1),....              17            0
#> 14         havana    core  Transcript c("ENSE0....              17            0
#> 15         havana    core  Transcript c("17", ....              17            0
#> 16         havana    core  Transcript c(1, 1, ....              17            0
#> 17         havana    core  Transcript c("ENSE0....              17            0
#> 18         havana    core  Transcript c("ENSE0....              17            0
#> 19         havana    core  Transcript c(-1, -1....              17            0
#> 20         havana    core  Transcript c(-1, -1....              17            0
#> 21         havana    core  Transcript c("core"....              17            0
#> 22         havana    core  Transcript c(412774....              17            0
#> 23         havana    core  Transcript c("core"....              17            0
#> 24         havana    core  Transcript c("Exon"....              17            0
#> 25         havana    core  Transcript c(412773....              17            0
#> 26         havana    core  Transcript c(412771....              17            0
#> 27        ensembl    core  Transcript c("GRCh3....              17            0
#> 28        ensembl    core  Transcript c("GRCh3....              17            0
#> 29        ensembl    core  Transcript c(412772....              17            0
#> 30        ensembl    core  Transcript c("ENSE0....              17            0
#> 31        ensembl    core  Transcript c("Exon"....              17            0
#>    species Translation.Parent  Translation.id Translation.db_type
#> 1    human    ENST00000357654 ENSP00000350283                core
#> 2    human    ENST00000468300 ENSP00000417148                core
#> 3    human    ENST00000586385 ENSP00000465818                core
#> 4    human    ENST00000591534 ENSP00000467329                core
#> 5    human    ENST00000591849 ENSP00000465347                core
#> 6    human    ENST00000493795 ENSP00000418775                core
#> 7    human    ENST00000471181 ENSP00000418960                core
#> 8    human    ENST00000461221 ENSP00000418548                core
#> 9    human    ENST00000491747 ENSP00000420705                core
#> 10   human    ENST00000484087 ENSP00000419481                core
#> 11   human    ENST00000478531 ENSP00000420412                core
#> 12   human    ENST00000493919 ENSP00000418819                core
#> 13   human               <NA>            <NA>                <NA>
#> 14   human    ENST00000487825 ENSP00000418212                core
#> 15   human    ENST00000461574 ENSP00000417241                core
#> 16   human               <NA>            <NA>                <NA>
#> 17   human    ENST00000412061 ENSP00000397145                core
#> 18   human    ENST00000470026 ENSP00000419274                core
#> 19   human    ENST00000477152 ENSP00000419988                core
#> 20   human    ENST00000492859 ENSP00000420253                core
#> 21   human    ENST00000497488 ENSP00000418986                core
#> 22   human    ENST00000494123 ENSP00000419103                core
#> 23   human    ENST00000473961 ENSP00000420201                core
#> 24   human    ENST00000476777 ENSP00000417554                core
#> 25   human    ENST00000461798 ENSP00000417988                core
#> 26   human    ENST00000489037 ENSP00000420781                core
#> 27   human    ENST00000354071 ENSP00000326002                core
#> 28   human    ENST00000352993 ENSP00000312236                core
#> 29   human    ENST00000346315 ENSP00000246907                core
#> 30   human    ENST00000351666 ENSP00000338007                core
#> 31   human    ENST00000309486 ENSP00000310938                core
#>    Translation.object_type Translation.start Translation.species
#> 1              Translation          41197695               human
#> 2              Translation          41197801               human
#> 3              Translation          41197695               human
#> 4              Translation          41197695               human
#> 5              Translation          41197695               human
#> 6              Translation          41197695               human
#> 7              Translation          41197695               human
#> 8              Translation          41256972               human
#> 9              Translation          41197695               human
#> 10             Translation          41215361               human
#> 11             Translation          41215361               human
#> 12             Translation          41215377               human
#> 13                    <NA>                NA                <NA>
#> 14             Translation          41228505               human
#> 15             Translation          41228554               human
#> 16                    <NA>                NA                <NA>
#> 17             Translation          41245587               human
#> 18             Translation          41245601               human
#> 19             Translation          41245603               human
#> 20             Translation          41262552               human
#> 21             Translation          41246129               human
#> 22             Translation          41246129               human
#> 23             Translation          41246187               human
#> 24             Translation          41247863               human
#> 25             Translation          41256972               human
#> 26             Translation          41256206               human
#> 27             Translation          41197695               human
#> 28             Translation          41197695               human
#> 29             Translation          41197695               human
#> 30             Translation          41197695               human
#> 31             Translation          41197695               human
#>    Translation.version Translation.length Translation.end strand
#> 1                    3               1863        41276113     -1
#> 2                    1                699        41276113     -1
#> 3                    1                173        41277202     -1
#> 4                    1                354        41226495     -1
#> 5                    1                 96        41202109     -1
#> 6                    1               1816        41258543     -1
#> 7                    2               1884        41276113     -1
#> 8                    1                 63        41276113     -1
#> 9                    2                759        41276113     -1
#> 10                   1                498        41256933     -1
#> 11                   1                623        41276113     -1
#> 12                   1                572        41258543     -1
#> 13                  NA                 NA              NA     -1
#> 14                   1                266        41256933     -1
#> 15                   1                242        41243841     -1
#> 16                  NA                 NA              NA     -1
#> 17                   3                437        41247883     -1
#> 18                   1                649        41276113     -1
#> 19                   1                622        41276113     -1
#> 20                   1                 59        41276113     -1
#> 21                   1                177        41246659     -1
#> 22                   1                473        41276113     -1
#> 23                   1                319        41256908     -1
#> 24                   1                222        41276113     -1
#> 25                   1                 63        41276113     -1
#> 26                   1                 98        41276113     -1
#> 27                   6               1598        41276113     -1
#> 28                   5                721        41276113     -1
#> 29                   4               1624        41276113     -1
#> 30                   3                680        41276113     -1
#> 31                   4               1567        41246659     -1
#>                    biotype                logic_name              id
#> 1           protein_coding ensembl_havana_transcript ENST00000357654
#> 2           protein_coding ensembl_havana_transcript ENST00000468300
#> 3           protein_coding    havana_homo_sapiens_37 ENST00000586385
#> 4           protein_coding    havana_homo_sapiens_37 ENST00000591534
#> 5           protein_coding    havana_homo_sapiens_37 ENST00000591849
#> 6           protein_coding ensembl_havana_transcript ENST00000493795
#> 7           protein_coding ensembl_havana_transcript ENST00000471181
#> 8  nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000461221
#> 9           protein_coding    havana_homo_sapiens_37 ENST00000491747
#> 10          protein_coding    havana_homo_sapiens_37 ENST00000484087
#> 11          protein_coding    havana_homo_sapiens_37 ENST00000478531
#> 12          protein_coding    havana_homo_sapiens_37 ENST00000493919
#> 13         retained_intron    havana_homo_sapiens_37 ENST00000472490
#> 14          protein_coding    havana_homo_sapiens_37 ENST00000487825
#> 15          protein_coding    havana_homo_sapiens_37 ENST00000461574
#> 16         retained_intron    havana_homo_sapiens_37 ENST00000467274
#> 17          non_stop_decay    havana_homo_sapiens_37 ENST00000412061
#> 18          protein_coding    havana_homo_sapiens_37 ENST00000470026
#> 19          protein_coding    havana_homo_sapiens_37 ENST00000477152
#> 20 nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000492859
#> 21          protein_coding    havana_homo_sapiens_37 ENST00000497488
#> 22          protein_coding    havana_homo_sapiens_37 ENST00000494123
#> 23          protein_coding    havana_homo_sapiens_37 ENST00000473961
#> 24          protein_coding    havana_homo_sapiens_37 ENST00000476777
#> 25 nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000461798
#> 26          protein_coding    havana_homo_sapiens_37 ENST00000489037
#> 27          protein_coding   ensembl_homo_sapiens_37 ENST00000354071
#> 28          protein_coding   ensembl_homo_sapiens_37 ENST00000352993
#> 29          protein_coding   ensembl_homo_sapiens_37 ENST00000346315
#> 30          protein_coding   ensembl_homo_sapiens_37 ENST00000351666
#> 31          protein_coding   ensembl_homo_sapiens_37 ENST00000309486
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
#> $species
#> [1] "human"
#> 
#> $seq_region_name
#> [1] "17"
#> 
```
