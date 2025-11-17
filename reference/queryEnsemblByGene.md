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
#> $db_type
#> [1] "core"
#> 
#> $end
#> [1] 41277500
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $species
#> [1] "human"
#> 
#> $strand
#> [1] -1
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
#> $version
#> [1] 15
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $Transcript
#>    strand         Exon species is_canonical      end seq_region_name
#> 1      -1 c("17", ....   human            0 41277387              17
#> 2      -1 c(412772....   human            0 41277468              17
#> 3      -1 c("17", ....   human            0 41277346              17
#> 4      -1 c(412772....   human            0 41277346              17
#> 5      -1 c(-1, -1....   human            0 41277346              17
#> 6      -1 c(412772....   human            0 41277419              17
#> 7      -1 c(412775....   human            1 41277500              17
#> 8      -1 c(-1, -1....   human            0 41277305              17
#> 9      -1 c(-1, -1....   human            0 41277373              17
#> 10     -1 c("ENSE0....   human            0 41256933              17
#> 11     -1 c("human....   human            0 41277376              17
#> 12     -1 c("17", ....   human            0 41277419              17
#> 13     -1 c("GRCh3....   human            0 41223083              17
#> 14     -1 c(-1, -1....   human            0 41256933              17
#> 15     -1 c("17", ....   human            0 41243841              17
#> 16     -1 c(412772....   human            0 41277332              17
#> 17     -1 c("17", ....   human            0 41247883              17
#> 18     -1 c("ENSE0....   human            0 41277340              17
#> 19     -1 c("human....   human            0 41277381              17
#> 20     -1 c("human....   human            0 41277317              17
#> 21     -1 c(412772....   human            0 41277317              17
#> 22     -1 c(-1, -1....   human            0 41277467              17
#> 23     -1 c("Exon"....   human            0 41256908              17
#> 24     -1 c("core"....   human            0 41277370              17
#> 25     -1 c(-1, -1....   human            0 41277387              17
#> 26     -1 c(412773....   human            0 41277338              17
#> 27     -1 c("ENSE0....   human            0 41277500              17
#> 28     -1 c(1, 1, ....   human            0 41277500              17
#> 29     -1 c(-1, -1....   human            0 41277468              17
#> 30     -1 c("17", ....   human            0 41276132              17
#> 31     -1 c(412772....   human            0 41277468              17
#>             Parent db_type         source length display_name
#> 1  ENSG00000012048    core ensembl_havana   7094    BRCA1-001
#> 2  ENSG00000012048    core ensembl_havana   3273    BRCA1-007
#> 3  ENSG00000012048    core         havana    781    BRCA1-023
#> 4  ENSG00000012048    core         havana   1282    BRCA1-024
#> 5  ENSG00000012048    core         havana    563    BRCA1-025
#> 6  ENSG00000012048    core ensembl_havana   5732    BRCA1-006
#> 7  ENSG00000012048    core ensembl_havana   5936    BRCA1-005
#> 8  ENSG00000012048    core         havana   5693    BRCA1-010
#> 9  ENSG00000012048    core         havana   2379    BRCA1-014
#> 10 ENSG00000012048    core         havana   1495    BRCA1-015
#> 11 ENSG00000012048    core         havana   1972    BRCA1-009
#> 12 ENSG00000012048    core         havana   1948    BRCA1-008
#> 13 ENSG00000012048    core         havana    561    BRCA1-021
#> 14 ENSG00000012048    core         havana    800    BRCA1-019
#> 15 ENSG00000012048    core         havana    726    BRCA1-022
#> 16 ENSG00000012048    core         havana   4497    BRCA1-012
#> 17 ENSG00000012048    core         havana   1312    BRCA1-026
#> 18 ENSG00000012048    core         havana   2108    BRCA1-011
#> 19 ENSG00000012048    core         havana   1980    BRCA1-004
#> 20 ENSG00000012048    core         havana   1584    BRCA1-002
#> 21 ENSG00000012048    core         havana    779    BRCA1-003
#> 22 ENSG00000012048    core         havana   1612    BRCA1-013
#> 23 ENSG00000012048    core         havana    958    BRCA1-018
#> 24 ENSG00000012048    core         havana    769    BRCA1-017
#> 25 ENSG00000012048    core         havana    582    BRCA1-020
#> 26 ENSG00000012048    core         havana    455    BRCA1-016
#> 27 ENSG00000012048    core        ensembl   6411    BRCA1-205
#> 28 ENSG00000012048    core        ensembl   3780    BRCA1-204
#> 29 ENSG00000012048    core        ensembl   6451    BRCA1-202
#> 30 ENSG00000012048    core        ensembl   3444    BRCA1-203
#> 31 ENSG00000012048    core        ensembl   7114    BRCA1-201
#>    Translation.length Translation.start Translation.species Translation.Parent
#> 1                1863          41197695               human    ENST00000357654
#> 2                 699          41197801               human    ENST00000468300
#> 3                 173          41197695               human    ENST00000586385
#> 4                 354          41197695               human    ENST00000591534
#> 5                  96          41197695               human    ENST00000591849
#> 6                1816          41197695               human    ENST00000493795
#> 7                1884          41197695               human    ENST00000471181
#> 8                  63          41256972               human    ENST00000461221
#> 9                 759          41197695               human    ENST00000491747
#> 10                498          41215361               human    ENST00000484087
#> 11                623          41215361               human    ENST00000478531
#> 12                572          41215377               human    ENST00000493919
#> 13                 NA                NA                <NA>               <NA>
#> 14                266          41228505               human    ENST00000487825
#> 15                242          41228554               human    ENST00000461574
#> 16                 NA                NA                <NA>               <NA>
#> 17                437          41245587               human    ENST00000412061
#> 18                649          41245601               human    ENST00000470026
#> 19                622          41245603               human    ENST00000477152
#> 20                 59          41262552               human    ENST00000492859
#> 21                177          41246129               human    ENST00000497488
#> 22                473          41246129               human    ENST00000494123
#> 23                319          41246187               human    ENST00000473961
#> 24                222          41247863               human    ENST00000476777
#> 25                 63          41256972               human    ENST00000461798
#> 26                 98          41256206               human    ENST00000489037
#> 27               1598          41197695               human    ENST00000354071
#> 28                721          41197695               human    ENST00000352993
#> 29               1624          41197695               human    ENST00000346315
#> 30                680          41197695               human    ENST00000351666
#> 31               1567          41197695               human    ENST00000309486
#>    Translation.end  Translation.id Translation.version Translation.db_type
#> 1         41276113 ENSP00000350283                   3                core
#> 2         41276113 ENSP00000417148                   1                core
#> 3         41277202 ENSP00000465818                   1                core
#> 4         41226495 ENSP00000467329                   1                core
#> 5         41202109 ENSP00000465347                   1                core
#> 6         41258543 ENSP00000418775                   1                core
#> 7         41276113 ENSP00000418960                   2                core
#> 8         41276113 ENSP00000418548                   1                core
#> 9         41276113 ENSP00000420705                   2                core
#> 10        41256933 ENSP00000419481                   1                core
#> 11        41276113 ENSP00000420412                   1                core
#> 12        41258543 ENSP00000418819                   1                core
#> 13              NA            <NA>                  NA                <NA>
#> 14        41256933 ENSP00000418212                   1                core
#> 15        41243841 ENSP00000417241                   1                core
#> 16              NA            <NA>                  NA                <NA>
#> 17        41247883 ENSP00000397145                   3                core
#> 18        41276113 ENSP00000419274                   1                core
#> 19        41276113 ENSP00000419988                   1                core
#> 20        41276113 ENSP00000420253                   1                core
#> 21        41246659 ENSP00000418986                   1                core
#> 22        41276113 ENSP00000419103                   1                core
#> 23        41256908 ENSP00000420201                   1                core
#> 24        41276113 ENSP00000417554                   1                core
#> 25        41276113 ENSP00000417988                   1                core
#> 26        41276113 ENSP00000420781                   1                core
#> 27        41276113 ENSP00000326002                   6                core
#> 28        41276113 ENSP00000312236                   5                core
#> 29        41276113 ENSP00000246907                   4                core
#> 30        41276113 ENSP00000338007                   3                core
#> 31        41246659 ENSP00000310938                   4                core
#>    Translation.object_type    start              id gencode_primary
#> 1              Translation 41196312 ENST00000357654               0
#> 2              Translation 41196822 ENST00000468300               0
#> 3              Translation 41197580 ENST00000586385               0
#> 4              Translation 41197580 ENST00000591534               0
#> 5              Translation 41197580 ENST00000591849               0
#> 6              Translation 41197646 ENST00000493795               0
#> 7              Translation 41197646 ENST00000471181               0
#> 8              Translation 41197695 ENST00000461221               0
#> 9              Translation 41197695 ENST00000491747               0
#> 10             Translation 41215361 ENST00000484087               0
#> 11             Translation 41215361 ENST00000478531               0
#> 12             Translation 41215377 ENST00000493919               0
#> 13                    <NA> 41219291 ENST00000472490               0
#> 14             Translation 41228505 ENST00000487825               0
#> 15             Translation 41228554 ENST00000461574               0
#> 16                    <NA> 41243115 ENST00000467274               0
#> 17             Translation 41245587 ENST00000412061               0
#> 18             Translation 41245601 ENST00000470026               0
#> 19             Translation 41245603 ENST00000477152               0
#> 20             Translation 41246129 ENST00000492859               0
#> 21             Translation 41246129 ENST00000497488               0
#> 22             Translation 41246129 ENST00000494123               0
#> 23             Translation 41246187 ENST00000473961               0
#> 24             Translation 41247863 ENST00000476777               0
#> 25             Translation 41251848 ENST00000461798               0
#> 26             Translation 41256206 ENST00000489037               0
#> 27             Translation 41196313 ENST00000354071               0
#> 28             Translation 41196313 ENST00000352993               0
#> 29             Translation 41196313 ENST00000346315               0
#> 30             Translation 41196313 ENST00000351666               0
#> 31             Translation 41196313 ENST00000309486               0
#>                   logic_name object_type                 biotype assembly_name
#> 1  ensembl_havana_transcript  Transcript          protein_coding        GRCh37
#> 2  ensembl_havana_transcript  Transcript          protein_coding        GRCh37
#> 3     havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 4     havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 5     havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 6  ensembl_havana_transcript  Transcript          protein_coding        GRCh37
#> 7  ensembl_havana_transcript  Transcript          protein_coding        GRCh37
#> 8     havana_homo_sapiens_37  Transcript nonsense_mediated_decay        GRCh37
#> 9     havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 10    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 11    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 12    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 13    havana_homo_sapiens_37  Transcript         retained_intron        GRCh37
#> 14    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 15    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 16    havana_homo_sapiens_37  Transcript         retained_intron        GRCh37
#> 17    havana_homo_sapiens_37  Transcript          non_stop_decay        GRCh37
#> 18    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 19    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 20    havana_homo_sapiens_37  Transcript nonsense_mediated_decay        GRCh37
#> 21    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 22    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 23    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 24    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 25    havana_homo_sapiens_37  Transcript nonsense_mediated_decay        GRCh37
#> 26    havana_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 27   ensembl_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 28   ensembl_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 29   ensembl_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 30   ensembl_homo_sapiens_37  Transcript          protein_coding        GRCh37
#> 31   ensembl_homo_sapiens_37  Transcript          protein_coding        GRCh37
#>    version
#> 1        3
#> 2        1
#> 3        1
#> 4        1
#> 5        1
#> 6        1
#> 7        2
#> 8        1
#> 9        2
#> 10       1
#> 11       1
#> 12       1
#> 13       1
#> 14       1
#> 15       1
#> 16       1
#> 17       3
#> 18       1
#> 19       1
#> 20       1
#> 21       1
#> 22       1
#> 23       1
#> 24       1
#> 25       1
#> 26       1
#> 27       3
#> 28       3
#> 29       3
#> 30       3
#> 31       4
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $start
#> [1] 41196312
#> 
#> $source
#> [1] "ensembl_havana"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $display_name
#> [1] "BRCA2"
#> 
#> $strand
#> [1] 1
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $start
#> [1] 32889611
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $end
#> [1] 32973805
#> 
#> $db_type
#> [1] "core"
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $version
#> [1] 10
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $Transcript
#>      start object_type gencode_primary strand display_name is_canonical
#> 1 32889611  Transcript               0      1    BRCA2-001            0
#> 2 32889642  Transcript               0      1    BRCA2-003            0
#> 3 32945108  Transcript               0      1    BRCA2-005            0
#> 4 32953977  Transcript               0      1    BRCA2-002            0
#> 5 32970946  Transcript               0      1    BRCA2-006            0
#> 6 32889617  Transcript               0      1    BRCA2-201            1
#>   Translation.object_type Translation.start Translation.Parent
#> 1             Translation          32890598    ENST00000380152
#> 2             Translation          32899266    ENST00000530893
#> 3             Translation          32945108    ENST00000528762
#> 4             Translation          32953977    ENST00000470094
#> 5                    <NA>                NA               <NA>
#> 6             Translation          32890598    ENST00000544455
#>   Translation.db_type Translation.end  Translation.id Translation.species
#> 1                core        32972907 ENSP00000369497        homo_sapiens
#> 2                core        32907428 ENSP00000435699        homo_sapiens
#> 3                core        32950807 ENSP00000433168        homo_sapiens
#> 4                core        32970229 ENSP00000434898        homo_sapiens
#> 5                <NA>              NA            <NA>                <NA>
#> 6                core        32972907 ENSP00000439902        homo_sapiens
#>   Translation.length Translation.version length db_type          Parent
#> 1               3418                   3  10930    core ENSG00000139618
#> 2                481                   2   2011    core ENSG00000139618
#> 3                 64                   1    495    core ENSG00000139618
#> 4                186                   1    842    core ENSG00000139618
#> 5                 NA                  NA    523    core ENSG00000139618
#> 6               3418                   1  10984    core ENSG00000139618
#>        end                logic_name         source              id
#> 1 32973347 ensembl_havana_transcript ensembl_havana ENST00000380152
#> 2 32907428    havana_homo_sapiens_37         havana ENST00000530893
#> 3 32953632    havana_homo_sapiens_37         havana ENST00000528762
#> 4 32972409    havana_homo_sapiens_37         havana ENST00000470094
#> 5 32972585    havana_homo_sapiens_37         havana ENST00000533776
#> 6 32973805   ensembl_homo_sapiens_37        ensembl ENST00000544455
#>        species         Exon version                 biotype assembly_name
#> 1 homo_sapiens c("13", ....       3          protein_coding        GRCh37
#> 2 homo_sapiens c("Exon"....       2          protein_coding        GRCh37
#> 3 homo_sapiens c("homo_....       1 nonsense_mediated_decay        GRCh37
#> 4 homo_sapiens c("ENSE0....       1 nonsense_mediated_decay        GRCh37
#> 5 homo_sapiens c("Exon"....       1         retained_intron        GRCh37
#> 6 homo_sapiens c(1, 1, ....       1          protein_coding        GRCh37
#>   seq_region_name
#> 1              13
#> 2              13
#> 3              13
#> 4              13
#> 5              13
#> 6              13
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $display_name
#> [1] "BRCA1"
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
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $version
#> [1] 15
#> 
#> $species
#> [1] "human"
#> 
#> $strand
#> [1] -1
#> 
#> $start
#> [1] 41196312
#> 
#> $end
#> [1] 41277500
#> 
#> $Transcript
#>    is_canonical species strand version display_name assembly_name
#> 1             0   human     -1       3    BRCA1-001        GRCh37
#> 2             0   human     -1       1    BRCA1-007        GRCh37
#> 3             0   human     -1       1    BRCA1-023        GRCh37
#> 4             0   human     -1       1    BRCA1-024        GRCh37
#> 5             0   human     -1       1    BRCA1-025        GRCh37
#> 6             0   human     -1       1    BRCA1-006        GRCh37
#> 7             1   human     -1       2    BRCA1-005        GRCh37
#> 8             0   human     -1       1    BRCA1-010        GRCh37
#> 9             0   human     -1       2    BRCA1-014        GRCh37
#> 10            0   human     -1       1    BRCA1-015        GRCh37
#> 11            0   human     -1       1    BRCA1-009        GRCh37
#> 12            0   human     -1       1    BRCA1-008        GRCh37
#> 13            0   human     -1       1    BRCA1-021        GRCh37
#> 14            0   human     -1       1    BRCA1-019        GRCh37
#> 15            0   human     -1       1    BRCA1-022        GRCh37
#> 16            0   human     -1       1    BRCA1-012        GRCh37
#> 17            0   human     -1       3    BRCA1-026        GRCh37
#> 18            0   human     -1       1    BRCA1-011        GRCh37
#> 19            0   human     -1       1    BRCA1-004        GRCh37
#> 20            0   human     -1       1    BRCA1-002        GRCh37
#> 21            0   human     -1       1    BRCA1-003        GRCh37
#> 22            0   human     -1       1    BRCA1-013        GRCh37
#> 23            0   human     -1       1    BRCA1-018        GRCh37
#> 24            0   human     -1       1    BRCA1-017        GRCh37
#> 25            0   human     -1       1    BRCA1-020        GRCh37
#> 26            0   human     -1       1    BRCA1-016        GRCh37
#> 27            0   human     -1       3    BRCA1-205        GRCh37
#> 28            0   human     -1       3    BRCA1-204        GRCh37
#> 29            0   human     -1       3    BRCA1-202        GRCh37
#> 30            0   human     -1       3    BRCA1-203        GRCh37
#> 31            0   human     -1       4    BRCA1-201        GRCh37
#>             Parent         source seq_region_name Translation.object_type
#> 1  ENSG00000012048 ensembl_havana              17             Translation
#> 2  ENSG00000012048 ensembl_havana              17             Translation
#> 3  ENSG00000012048         havana              17             Translation
#> 4  ENSG00000012048         havana              17             Translation
#> 5  ENSG00000012048         havana              17             Translation
#> 6  ENSG00000012048 ensembl_havana              17             Translation
#> 7  ENSG00000012048 ensembl_havana              17             Translation
#> 8  ENSG00000012048         havana              17             Translation
#> 9  ENSG00000012048         havana              17             Translation
#> 10 ENSG00000012048         havana              17             Translation
#> 11 ENSG00000012048         havana              17             Translation
#> 12 ENSG00000012048         havana              17             Translation
#> 13 ENSG00000012048         havana              17                    <NA>
#> 14 ENSG00000012048         havana              17             Translation
#> 15 ENSG00000012048         havana              17             Translation
#> 16 ENSG00000012048         havana              17                    <NA>
#> 17 ENSG00000012048         havana              17             Translation
#> 18 ENSG00000012048         havana              17             Translation
#> 19 ENSG00000012048         havana              17             Translation
#> 20 ENSG00000012048         havana              17             Translation
#> 21 ENSG00000012048         havana              17             Translation
#> 22 ENSG00000012048         havana              17             Translation
#> 23 ENSG00000012048         havana              17             Translation
#> 24 ENSG00000012048         havana              17             Translation
#> 25 ENSG00000012048         havana              17             Translation
#> 26 ENSG00000012048         havana              17             Translation
#> 27 ENSG00000012048        ensembl              17             Translation
#> 28 ENSG00000012048        ensembl              17             Translation
#> 29 ENSG00000012048        ensembl              17             Translation
#> 30 ENSG00000012048        ensembl              17             Translation
#> 31 ENSG00000012048        ensembl              17             Translation
#>    Translation.Parent  Translation.id Translation.length Translation.db_type
#> 1     ENST00000357654 ENSP00000350283               1863                core
#> 2     ENST00000468300 ENSP00000417148                699                core
#> 3     ENST00000586385 ENSP00000465818                173                core
#> 4     ENST00000591534 ENSP00000467329                354                core
#> 5     ENST00000591849 ENSP00000465347                 96                core
#> 6     ENST00000493795 ENSP00000418775               1816                core
#> 7     ENST00000471181 ENSP00000418960               1884                core
#> 8     ENST00000461221 ENSP00000418548                 63                core
#> 9     ENST00000491747 ENSP00000420705                759                core
#> 10    ENST00000484087 ENSP00000419481                498                core
#> 11    ENST00000478531 ENSP00000420412                623                core
#> 12    ENST00000493919 ENSP00000418819                572                core
#> 13               <NA>            <NA>                 NA                <NA>
#> 14    ENST00000487825 ENSP00000418212                266                core
#> 15    ENST00000461574 ENSP00000417241                242                core
#> 16               <NA>            <NA>                 NA                <NA>
#> 17    ENST00000412061 ENSP00000397145                437                core
#> 18    ENST00000470026 ENSP00000419274                649                core
#> 19    ENST00000477152 ENSP00000419988                622                core
#> 20    ENST00000492859 ENSP00000420253                 59                core
#> 21    ENST00000497488 ENSP00000418986                177                core
#> 22    ENST00000494123 ENSP00000419103                473                core
#> 23    ENST00000473961 ENSP00000420201                319                core
#> 24    ENST00000476777 ENSP00000417554                222                core
#> 25    ENST00000461798 ENSP00000417988                 63                core
#> 26    ENST00000489037 ENSP00000420781                 98                core
#> 27    ENST00000354071 ENSP00000326002               1598                core
#> 28    ENST00000352993 ENSP00000312236                721                core
#> 29    ENST00000346315 ENSP00000246907               1624                core
#> 30    ENST00000351666 ENSP00000338007                680                core
#> 31    ENST00000309486 ENSP00000310938               1567                core
#>    Translation.end Translation.start Translation.species Translation.version
#> 1         41276113          41197695               human                   3
#> 2         41276113          41197801               human                   1
#> 3         41277202          41197695               human                   1
#> 4         41226495          41197695               human                   1
#> 5         41202109          41197695               human                   1
#> 6         41258543          41197695               human                   1
#> 7         41276113          41197695               human                   2
#> 8         41276113          41256972               human                   1
#> 9         41276113          41197695               human                   2
#> 10        41256933          41215361               human                   1
#> 11        41276113          41215361               human                   1
#> 12        41258543          41215377               human                   1
#> 13              NA                NA                <NA>                  NA
#> 14        41256933          41228505               human                   1
#> 15        41243841          41228554               human                   1
#> 16              NA                NA                <NA>                  NA
#> 17        41247883          41245587               human                   3
#> 18        41276113          41245601               human                   1
#> 19        41276113          41245603               human                   1
#> 20        41276113          41262552               human                   1
#> 21        41246659          41246129               human                   1
#> 22        41276113          41246129               human                   1
#> 23        41256908          41246187               human                   1
#> 24        41276113          41247863               human                   1
#> 25        41276113          41256972               human                   1
#> 26        41276113          41256206               human                   1
#> 27        41276113          41197695               human                   6
#> 28        41276113          41197695               human                   5
#> 29        41276113          41197695               human                   4
#> 30        41276113          41197695               human                   3
#> 31        41246659          41197695               human                   4
#>                 id                logic_name gencode_primary         Exon
#> 1  ENST00000357654 ensembl_havana_transcript               0 c("GRCh3....
#> 2  ENST00000468300 ensembl_havana_transcript               0 c(412774....
#> 3  ENST00000586385    havana_homo_sapiens_37               0 c(1, 1, ....
#> 4  ENST00000591534    havana_homo_sapiens_37               0 c("core"....
#> 5  ENST00000591849    havana_homo_sapiens_37               0 c(412773....
#> 6  ENST00000493795 ensembl_havana_transcript               0 c("Exon"....
#> 7  ENST00000471181 ensembl_havana_transcript               0 c(412772....
#> 8  ENST00000461221    havana_homo_sapiens_37               0 c("core"....
#> 9  ENST00000491747    havana_homo_sapiens_37               0 c(412772....
#> 10 ENST00000484087    havana_homo_sapiens_37               0 c("Exon"....
#> 11 ENST00000478531    havana_homo_sapiens_37               0 c(1, 1, ....
#> 12 ENST00000493919    havana_homo_sapiens_37               0 c(412772....
#> 13 ENST00000472490    havana_homo_sapiens_37               0 c("17", ....
#> 14 ENST00000487825    havana_homo_sapiens_37               0 c("17", ....
#> 15 ENST00000461574    havana_homo_sapiens_37               0 c("human....
#> 16 ENST00000467274    havana_homo_sapiens_37               0 c("human....
#> 17 ENST00000412061    havana_homo_sapiens_37               0 c(412478....
#> 18 ENST00000470026    havana_homo_sapiens_37               0 c("human....
#> 19 ENST00000477152    havana_homo_sapiens_37               0 c("ENSE0....
#> 20 ENST00000492859    havana_homo_sapiens_37               0 c("Exon"....
#> 21 ENST00000497488    havana_homo_sapiens_37               0 c("human....
#> 22 ENST00000494123    havana_homo_sapiens_37               0 c("ENSE0....
#> 23 ENST00000473961    havana_homo_sapiens_37               0 c("Exon"....
#> 24 ENST00000476777    havana_homo_sapiens_37               0 c(412773....
#> 25 ENST00000461798    havana_homo_sapiens_37               0 c("human....
#> 26 ENST00000489037    havana_homo_sapiens_37               0 c(412771....
#> 27 ENST00000354071   ensembl_homo_sapiens_37               0 c("ENSE0....
#> 28 ENST00000352993   ensembl_homo_sapiens_37               0 c("Exon"....
#> 29 ENST00000346315   ensembl_homo_sapiens_37               0 c(-1, -1....
#> 30 ENST00000351666   ensembl_homo_sapiens_37               0 c(-1, -1....
#> 31 ENST00000309486   ensembl_homo_sapiens_37               0 c(412774....
#>         end db_type    start object_type                 biotype length
#> 1  41277387    core 41196312  Transcript          protein_coding   7094
#> 2  41277468    core 41196822  Transcript          protein_coding   3273
#> 3  41277346    core 41197580  Transcript          protein_coding    781
#> 4  41277346    core 41197580  Transcript          protein_coding   1282
#> 5  41277346    core 41197580  Transcript          protein_coding    563
#> 6  41277419    core 41197646  Transcript          protein_coding   5732
#> 7  41277500    core 41197646  Transcript          protein_coding   5936
#> 8  41277305    core 41197695  Transcript nonsense_mediated_decay   5693
#> 9  41277373    core 41197695  Transcript          protein_coding   2379
#> 10 41256933    core 41215361  Transcript          protein_coding   1495
#> 11 41277376    core 41215361  Transcript          protein_coding   1972
#> 12 41277419    core 41215377  Transcript          protein_coding   1948
#> 13 41223083    core 41219291  Transcript         retained_intron    561
#> 14 41256933    core 41228505  Transcript          protein_coding    800
#> 15 41243841    core 41228554  Transcript          protein_coding    726
#> 16 41277332    core 41243115  Transcript         retained_intron   4497
#> 17 41247883    core 41245587  Transcript          non_stop_decay   1312
#> 18 41277340    core 41245601  Transcript          protein_coding   2108
#> 19 41277381    core 41245603  Transcript          protein_coding   1980
#> 20 41277317    core 41246129  Transcript nonsense_mediated_decay   1584
#> 21 41277317    core 41246129  Transcript          protein_coding    779
#> 22 41277467    core 41246129  Transcript          protein_coding   1612
#> 23 41256908    core 41246187  Transcript          protein_coding    958
#> 24 41277370    core 41247863  Transcript          protein_coding    769
#> 25 41277387    core 41251848  Transcript nonsense_mediated_decay    582
#> 26 41277338    core 41256206  Transcript          protein_coding    455
#> 27 41277500    core 41196313  Transcript          protein_coding   6411
#> 28 41277500    core 41196313  Transcript          protein_coding   3780
#> 29 41277468    core 41196313  Transcript          protein_coding   6451
#> 30 41276132    core 41196313  Transcript          protein_coding   3444
#> 31 41277468    core 41196313  Transcript          protein_coding   7114
#> 
#> $db_type
#> [1] "core"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
```
