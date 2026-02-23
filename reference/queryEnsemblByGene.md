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
#> $display_name
#> [1] "BRCA1"
#> 
#> $start
#> [1] 41196312
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
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
#> $end
#> [1] 41277500
#> 
#> $db_type
#> [1] "core"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
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
#> $seq_region_name
#> [1] "17"
#> 
#> $Transcript
#>    db_type          Parent      end              id                logic_name
#> 1     core ENSG00000012048 41277387 ENST00000357654 ensembl_havana_transcript
#> 2     core ENSG00000012048 41277468 ENST00000468300 ensembl_havana_transcript
#> 3     core ENSG00000012048 41277346 ENST00000586385    havana_homo_sapiens_37
#> 4     core ENSG00000012048 41277346 ENST00000591534    havana_homo_sapiens_37
#> 5     core ENSG00000012048 41277346 ENST00000591849    havana_homo_sapiens_37
#> 6     core ENSG00000012048 41277419 ENST00000493795 ensembl_havana_transcript
#> 7     core ENSG00000012048 41277500 ENST00000471181 ensembl_havana_transcript
#> 8     core ENSG00000012048 41277305 ENST00000461221    havana_homo_sapiens_37
#> 9     core ENSG00000012048 41277373 ENST00000491747    havana_homo_sapiens_37
#> 10    core ENSG00000012048 41256933 ENST00000484087    havana_homo_sapiens_37
#> 11    core ENSG00000012048 41277376 ENST00000478531    havana_homo_sapiens_37
#> 12    core ENSG00000012048 41277419 ENST00000493919    havana_homo_sapiens_37
#> 13    core ENSG00000012048 41223083 ENST00000472490    havana_homo_sapiens_37
#> 14    core ENSG00000012048 41256933 ENST00000487825    havana_homo_sapiens_37
#> 15    core ENSG00000012048 41243841 ENST00000461574    havana_homo_sapiens_37
#> 16    core ENSG00000012048 41277332 ENST00000467274    havana_homo_sapiens_37
#> 17    core ENSG00000012048 41247883 ENST00000412061    havana_homo_sapiens_37
#> 18    core ENSG00000012048 41277340 ENST00000470026    havana_homo_sapiens_37
#> 19    core ENSG00000012048 41277381 ENST00000477152    havana_homo_sapiens_37
#> 20    core ENSG00000012048 41277317 ENST00000492859    havana_homo_sapiens_37
#> 21    core ENSG00000012048 41277317 ENST00000497488    havana_homo_sapiens_37
#> 22    core ENSG00000012048 41277467 ENST00000494123    havana_homo_sapiens_37
#> 23    core ENSG00000012048 41256908 ENST00000473961    havana_homo_sapiens_37
#> 24    core ENSG00000012048 41277370 ENST00000476777    havana_homo_sapiens_37
#> 25    core ENSG00000012048 41277387 ENST00000461798    havana_homo_sapiens_37
#> 26    core ENSG00000012048 41277338 ENST00000489037    havana_homo_sapiens_37
#> 27    core ENSG00000012048 41277500 ENST00000354071   ensembl_homo_sapiens_37
#> 28    core ENSG00000012048 41277500 ENST00000352993   ensembl_homo_sapiens_37
#> 29    core ENSG00000012048 41277468 ENST00000346315   ensembl_homo_sapiens_37
#> 30    core ENSG00000012048 41276132 ENST00000351666   ensembl_homo_sapiens_37
#> 31    core ENSG00000012048 41277468 ENST00000309486   ensembl_homo_sapiens_37
#>            source         Exon species seq_region_name                 biotype
#> 1  ensembl_havana c("Exon"....   human              17          protein_coding
#> 2  ensembl_havana c(-1, -1....   human              17          protein_coding
#> 3          havana c("ENSE0....   human              17          protein_coding
#> 4          havana c("17", ....   human              17          protein_coding
#> 5          havana c("human....   human              17          protein_coding
#> 6  ensembl_havana c("ENSE0....   human              17          protein_coding
#> 7  ensembl_havana c(-1, -1....   human              17          protein_coding
#> 8          havana c("17", ....   human              17 nonsense_mediated_decay
#> 9          havana c("ENSE0....   human              17          protein_coding
#> 10         havana c("Exon"....   human              17          protein_coding
#> 11         havana c(1, 1, ....   human              17          protein_coding
#> 12         havana c(-1, -1....   human              17          protein_coding
#> 13         havana c(1, 1),....   human              17         retained_intron
#> 14         havana c("ENSE0....   human              17          protein_coding
#> 15         havana c("GRCh3....   human              17          protein_coding
#> 16         havana c("human....   human              17         retained_intron
#> 17         havana c("core"....   human              17          non_stop_decay
#> 18         havana c("17", ....   human              17          protein_coding
#> 19         havana c(412772....   human              17          protein_coding
#> 20         havana c("ENSE0....   human              17 nonsense_mediated_decay
#> 21         havana c(-1, -1....   human              17          protein_coding
#> 22         havana c("human....   human              17          protein_coding
#> 23         havana c("human....   human              17          protein_coding
#> 24         havana c("core"....   human              17          protein_coding
#> 25         havana c("human....   human              17 nonsense_mediated_decay
#> 26         havana c("ENSE0....   human              17          protein_coding
#> 27        ensembl c(-1, -1....   human              17          protein_coding
#> 28        ensembl c(-1, -1....   human              17          protein_coding
#> 29        ensembl c("ENSE0....   human              17          protein_coding
#> 30        ensembl c("human....   human              17          protein_coding
#> 31        ensembl c("human....   human              17          protein_coding
#>    assembly_name version gencode_primary object_type    start display_name
#> 1         GRCh37       3               0  Transcript 41196312    BRCA1-001
#> 2         GRCh37       1               0  Transcript 41196822    BRCA1-007
#> 3         GRCh37       1               0  Transcript 41197580    BRCA1-023
#> 4         GRCh37       1               0  Transcript 41197580    BRCA1-024
#> 5         GRCh37       1               0  Transcript 41197580    BRCA1-025
#> 6         GRCh37       1               0  Transcript 41197646    BRCA1-006
#> 7         GRCh37       2               0  Transcript 41197646    BRCA1-005
#> 8         GRCh37       1               0  Transcript 41197695    BRCA1-010
#> 9         GRCh37       2               0  Transcript 41197695    BRCA1-014
#> 10        GRCh37       1               0  Transcript 41215361    BRCA1-015
#> 11        GRCh37       1               0  Transcript 41215361    BRCA1-009
#> 12        GRCh37       1               0  Transcript 41215377    BRCA1-008
#> 13        GRCh37       1               0  Transcript 41219291    BRCA1-021
#> 14        GRCh37       1               0  Transcript 41228505    BRCA1-019
#> 15        GRCh37       1               0  Transcript 41228554    BRCA1-022
#> 16        GRCh37       1               0  Transcript 41243115    BRCA1-012
#> 17        GRCh37       3               0  Transcript 41245587    BRCA1-026
#> 18        GRCh37       1               0  Transcript 41245601    BRCA1-011
#> 19        GRCh37       1               0  Transcript 41245603    BRCA1-004
#> 20        GRCh37       1               0  Transcript 41246129    BRCA1-002
#> 21        GRCh37       1               0  Transcript 41246129    BRCA1-003
#> 22        GRCh37       1               0  Transcript 41246129    BRCA1-013
#> 23        GRCh37       1               0  Transcript 41246187    BRCA1-018
#> 24        GRCh37       1               0  Transcript 41247863    BRCA1-017
#> 25        GRCh37       1               0  Transcript 41251848    BRCA1-020
#> 26        GRCh37       1               0  Transcript 41256206    BRCA1-016
#> 27        GRCh37       3               0  Transcript 41196313    BRCA1-205
#> 28        GRCh37       3               0  Transcript 41196313    BRCA1-204
#> 29        GRCh37       3               0  Transcript 41196313    BRCA1-202
#> 30        GRCh37       3               0  Transcript 41196313    BRCA1-203
#> 31        GRCh37       4               0  Transcript 41196313    BRCA1-201
#>    strand is_canonical length Translation.species Translation.length
#> 1      -1            0   7094               human               1863
#> 2      -1            0   3273               human                699
#> 3      -1            0    781               human                173
#> 4      -1            0   1282               human                354
#> 5      -1            0    563               human                 96
#> 6      -1            0   5732               human               1816
#> 7      -1            1   5936               human               1884
#> 8      -1            0   5693               human                 63
#> 9      -1            0   2379               human                759
#> 10     -1            0   1495               human                498
#> 11     -1            0   1972               human                623
#> 12     -1            0   1948               human                572
#> 13     -1            0    561                <NA>                 NA
#> 14     -1            0    800               human                266
#> 15     -1            0    726               human                242
#> 16     -1            0   4497                <NA>                 NA
#> 17     -1            0   1312               human                437
#> 18     -1            0   2108               human                649
#> 19     -1            0   1980               human                622
#> 20     -1            0   1584               human                 59
#> 21     -1            0    779               human                177
#> 22     -1            0   1612               human                473
#> 23     -1            0    958               human                319
#> 24     -1            0    769               human                222
#> 25     -1            0    582               human                 63
#> 26     -1            0    455               human                 98
#> 27     -1            0   6411               human               1598
#> 28     -1            0   3780               human                721
#> 29     -1            0   6451               human               1624
#> 30     -1            0   3444               human                680
#> 31     -1            0   7114               human               1567
#>    Translation.version Translation.object_type Translation.Parent
#> 1                    3             Translation    ENST00000357654
#> 2                    1             Translation    ENST00000468300
#> 3                    1             Translation    ENST00000586385
#> 4                    1             Translation    ENST00000591534
#> 5                    1             Translation    ENST00000591849
#> 6                    1             Translation    ENST00000493795
#> 7                    2             Translation    ENST00000471181
#> 8                    1             Translation    ENST00000461221
#> 9                    2             Translation    ENST00000491747
#> 10                   1             Translation    ENST00000484087
#> 11                   1             Translation    ENST00000478531
#> 12                   1             Translation    ENST00000493919
#> 13                  NA                    <NA>               <NA>
#> 14                   1             Translation    ENST00000487825
#> 15                   1             Translation    ENST00000461574
#> 16                  NA                    <NA>               <NA>
#> 17                   3             Translation    ENST00000412061
#> 18                   1             Translation    ENST00000470026
#> 19                   1             Translation    ENST00000477152
#> 20                   1             Translation    ENST00000492859
#> 21                   1             Translation    ENST00000497488
#> 22                   1             Translation    ENST00000494123
#> 23                   1             Translation    ENST00000473961
#> 24                   1             Translation    ENST00000476777
#> 25                   1             Translation    ENST00000461798
#> 26                   1             Translation    ENST00000489037
#> 27                   6             Translation    ENST00000354071
#> 28                   5             Translation    ENST00000352993
#> 29                   4             Translation    ENST00000346315
#> 30                   3             Translation    ENST00000351666
#> 31                   4             Translation    ENST00000309486
#>    Translation.db_type Translation.start Translation.end  Translation.id
#> 1                 core          41197695        41276113 ENSP00000350283
#> 2                 core          41197801        41276113 ENSP00000417148
#> 3                 core          41197695        41277202 ENSP00000465818
#> 4                 core          41197695        41226495 ENSP00000467329
#> 5                 core          41197695        41202109 ENSP00000465347
#> 6                 core          41197695        41258543 ENSP00000418775
#> 7                 core          41197695        41276113 ENSP00000418960
#> 8                 core          41256972        41276113 ENSP00000418548
#> 9                 core          41197695        41276113 ENSP00000420705
#> 10                core          41215361        41256933 ENSP00000419481
#> 11                core          41215361        41276113 ENSP00000420412
#> 12                core          41215377        41258543 ENSP00000418819
#> 13                <NA>                NA              NA            <NA>
#> 14                core          41228505        41256933 ENSP00000418212
#> 15                core          41228554        41243841 ENSP00000417241
#> 16                <NA>                NA              NA            <NA>
#> 17                core          41245587        41247883 ENSP00000397145
#> 18                core          41245601        41276113 ENSP00000419274
#> 19                core          41245603        41276113 ENSP00000419988
#> 20                core          41262552        41276113 ENSP00000420253
#> 21                core          41246129        41246659 ENSP00000418986
#> 22                core          41246129        41276113 ENSP00000419103
#> 23                core          41246187        41256908 ENSP00000420201
#> 24                core          41247863        41276113 ENSP00000417554
#> 25                core          41256972        41276113 ENSP00000417988
#> 26                core          41256206        41276113 ENSP00000420781
#> 27                core          41197695        41276113 ENSP00000326002
#> 28                core          41197695        41276113 ENSP00000312236
#> 29                core          41197695        41276113 ENSP00000246907
#> 30                core          41197695        41276113 ENSP00000338007
#> 31                core          41197695        41246659 ENSP00000310938
#> 
#> $species
#> [1] "human"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $start
#> [1] 32889611
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $strand
#> [1] 1
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $version
#> [1] 10
#> 
#> $db_type
#> [1] "core"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $Transcript
#>        end object_type display_name      species         Exon assembly_name
#> 1 32973347  Transcript    BRCA2-001 homo_sapiens c(328896....        GRCh37
#> 2 32907428  Transcript    BRCA2-003 homo_sapiens c(328898....        GRCh37
#> 3 32953632  Transcript    BRCA2-005 homo_sapiens c(329451....        GRCh37
#> 4 32972409  Transcript    BRCA2-002 homo_sapiens c("Exon"....        GRCh37
#> 5 32972585  Transcript    BRCA2-006 homo_sapiens c("core"....        GRCh37
#> 6 32973805  Transcript    BRCA2-201 homo_sapiens c("GRCh3....        GRCh37
#>   db_type is_canonical                logic_name              id strand
#> 1    core            0 ensembl_havana_transcript ENST00000380152      1
#> 2    core            0    havana_homo_sapiens_37 ENST00000530893      1
#> 3    core            0    havana_homo_sapiens_37 ENST00000528762      1
#> 4    core            0    havana_homo_sapiens_37 ENST00000470094      1
#> 5    core            0    havana_homo_sapiens_37 ENST00000533776      1
#> 6    core            1   ensembl_homo_sapiens_37 ENST00000544455      1
#>   seq_region_name version    start         source length gencode_primary
#> 1              13       3 32889611 ensembl_havana  10930               0
#> 2              13       2 32889642         havana   2011               0
#> 3              13       1 32945108         havana    495               0
#> 4              13       1 32953977         havana    842               0
#> 5              13       1 32970946         havana    523               0
#> 6              13       1 32889617        ensembl  10984               0
#>            Parent  Translation.id Translation.version Translation.db_type
#> 1 ENSG00000139618 ENSP00000369497                   3                core
#> 2 ENSG00000139618 ENSP00000435699                   2                core
#> 3 ENSG00000139618 ENSP00000433168                   1                core
#> 4 ENSG00000139618 ENSP00000434898                   1                core
#> 5 ENSG00000139618            <NA>                  NA                <NA>
#> 6 ENSG00000139618 ENSP00000439902                   1                core
#>   Translation.Parent Translation.start Translation.species
#> 1    ENST00000380152          32890598        homo_sapiens
#> 2    ENST00000530893          32899266        homo_sapiens
#> 3    ENST00000528762          32945108        homo_sapiens
#> 4    ENST00000470094          32953977        homo_sapiens
#> 5               <NA>                NA                <NA>
#> 6    ENST00000544455          32890598        homo_sapiens
#>   Translation.object_type Translation.end Translation.length
#> 1             Translation        32972907               3418
#> 2             Translation        32907428                481
#> 3             Translation        32950807                 64
#> 4             Translation        32970229                186
#> 5                    <NA>              NA                 NA
#> 6             Translation        32972907               3418
#>                   biotype
#> 1          protein_coding
#> 2          protein_coding
#> 3 nonsense_mediated_decay
#> 4 nonsense_mediated_decay
#> 5         retained_intron
#> 6          protein_coding
#> 
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $end
#> [1] 32973805
#> 
#> $display_name
#> [1] "BRCA2"
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $species
#> [1] "human"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $end
#> [1] 41277500
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $Transcript
#>    is_canonical db_type                logic_name strand version
#> 1             0    core ensembl_havana_transcript     -1       3
#> 2             0    core ensembl_havana_transcript     -1       1
#> 3             0    core    havana_homo_sapiens_37     -1       1
#> 4             0    core    havana_homo_sapiens_37     -1       1
#> 5             0    core    havana_homo_sapiens_37     -1       1
#> 6             0    core ensembl_havana_transcript     -1       1
#> 7             1    core ensembl_havana_transcript     -1       2
#> 8             0    core    havana_homo_sapiens_37     -1       1
#> 9             0    core    havana_homo_sapiens_37     -1       2
#> 10            0    core    havana_homo_sapiens_37     -1       1
#> 11            0    core    havana_homo_sapiens_37     -1       1
#> 12            0    core    havana_homo_sapiens_37     -1       1
#> 13            0    core    havana_homo_sapiens_37     -1       1
#> 14            0    core    havana_homo_sapiens_37     -1       1
#> 15            0    core    havana_homo_sapiens_37     -1       1
#> 16            0    core    havana_homo_sapiens_37     -1       1
#> 17            0    core    havana_homo_sapiens_37     -1       3
#> 18            0    core    havana_homo_sapiens_37     -1       1
#> 19            0    core    havana_homo_sapiens_37     -1       1
#> 20            0    core    havana_homo_sapiens_37     -1       1
#> 21            0    core    havana_homo_sapiens_37     -1       1
#> 22            0    core    havana_homo_sapiens_37     -1       1
#> 23            0    core    havana_homo_sapiens_37     -1       1
#> 24            0    core    havana_homo_sapiens_37     -1       1
#> 25            0    core    havana_homo_sapiens_37     -1       1
#> 26            0    core    havana_homo_sapiens_37     -1       1
#> 27            0    core   ensembl_homo_sapiens_37     -1       3
#> 28            0    core   ensembl_homo_sapiens_37     -1       3
#> 29            0    core   ensembl_homo_sapiens_37     -1       3
#> 30            0    core   ensembl_homo_sapiens_37     -1       3
#> 31            0    core   ensembl_homo_sapiens_37     -1       4
#>    seq_region_name              id assembly_name         Exon display_name
#> 1               17 ENST00000357654        GRCh37 c("core"....    BRCA1-001
#> 2               17 ENST00000468300        GRCh37 c("human....    BRCA1-007
#> 3               17 ENST00000586385        GRCh37 c("human....    BRCA1-023
#> 4               17 ENST00000591534        GRCh37 c("Exon"....    BRCA1-024
#> 5               17 ENST00000591849        GRCh37 c("Exon"....    BRCA1-025
#> 6               17 ENST00000493795        GRCh37 c("GRCh3....    BRCA1-006
#> 7               17 ENST00000471181        GRCh37 c("GRCh3....    BRCA1-005
#> 8               17 ENST00000461221        GRCh37 c(412773....    BRCA1-010
#> 9               17 ENST00000491747        GRCh37 c("core"....    BRCA1-014
#> 10              17 ENST00000484087        GRCh37 c("Exon"....    BRCA1-015
#> 11              17 ENST00000478531        GRCh37 c("human....    BRCA1-009
#> 12              17 ENST00000493919        GRCh37 c("human....    BRCA1-008
#> 13              17 ENST00000472490        GRCh37 c(412229....    BRCA1-021
#> 14              17 ENST00000487825        GRCh37 c("human....    BRCA1-019
#> 15              17 ENST00000461574        GRCh37 c("ENSE0....    BRCA1-022
#> 16              17 ENST00000467274        GRCh37 c("17", ....    BRCA1-012
#> 17              17 ENST00000412061        GRCh37 c("Exon"....    BRCA1-026
#> 18              17 ENST00000470026        GRCh37 c("human....    BRCA1-011
#> 19              17 ENST00000477152        GRCh37 c("core"....    BRCA1-004
#> 20              17 ENST00000492859        GRCh37 c("Exon"....    BRCA1-002
#> 21              17 ENST00000497488        GRCh37 c("GRCh3....    BRCA1-003
#> 22              17 ENST00000494123        GRCh37 c(412772....    BRCA1-013
#> 23              17 ENST00000473961        GRCh37 c("GRCh3....    BRCA1-018
#> 24              17 ENST00000476777        GRCh37 c("human....    BRCA1-017
#> 25              17 ENST00000461798        GRCh37 c("ENSE0....    BRCA1-020
#> 26              17 ENST00000489037        GRCh37 c(412771....    BRCA1-016
#> 27              17 ENST00000354071        GRCh37 c("human....    BRCA1-205
#> 28              17 ENST00000352993        GRCh37 c("GRCh3....    BRCA1-204
#> 29              17 ENST00000346315        GRCh37 c(412772....    BRCA1-202
#> 30              17 ENST00000351666        GRCh37 c("human....    BRCA1-203
#> 31              17 ENST00000309486        GRCh37 c("GRCh3....    BRCA1-201
#>    object_type      end species          Parent gencode_primary
#> 1   Transcript 41277387   human ENSG00000012048               0
#> 2   Transcript 41277468   human ENSG00000012048               0
#> 3   Transcript 41277346   human ENSG00000012048               0
#> 4   Transcript 41277346   human ENSG00000012048               0
#> 5   Transcript 41277346   human ENSG00000012048               0
#> 6   Transcript 41277419   human ENSG00000012048               0
#> 7   Transcript 41277500   human ENSG00000012048               0
#> 8   Transcript 41277305   human ENSG00000012048               0
#> 9   Transcript 41277373   human ENSG00000012048               0
#> 10  Transcript 41256933   human ENSG00000012048               0
#> 11  Transcript 41277376   human ENSG00000012048               0
#> 12  Transcript 41277419   human ENSG00000012048               0
#> 13  Transcript 41223083   human ENSG00000012048               0
#> 14  Transcript 41256933   human ENSG00000012048               0
#> 15  Transcript 41243841   human ENSG00000012048               0
#> 16  Transcript 41277332   human ENSG00000012048               0
#> 17  Transcript 41247883   human ENSG00000012048               0
#> 18  Transcript 41277340   human ENSG00000012048               0
#> 19  Transcript 41277381   human ENSG00000012048               0
#> 20  Transcript 41277317   human ENSG00000012048               0
#> 21  Transcript 41277317   human ENSG00000012048               0
#> 22  Transcript 41277467   human ENSG00000012048               0
#> 23  Transcript 41256908   human ENSG00000012048               0
#> 24  Transcript 41277370   human ENSG00000012048               0
#> 25  Transcript 41277387   human ENSG00000012048               0
#> 26  Transcript 41277338   human ENSG00000012048               0
#> 27  Transcript 41277500   human ENSG00000012048               0
#> 28  Transcript 41277500   human ENSG00000012048               0
#> 29  Transcript 41277468   human ENSG00000012048               0
#> 30  Transcript 41276132   human ENSG00000012048               0
#> 31  Transcript 41277468   human ENSG00000012048               0
#>                    biotype Translation.start Translation.species
#> 1           protein_coding          41197695               human
#> 2           protein_coding          41197801               human
#> 3           protein_coding          41197695               human
#> 4           protein_coding          41197695               human
#> 5           protein_coding          41197695               human
#> 6           protein_coding          41197695               human
#> 7           protein_coding          41197695               human
#> 8  nonsense_mediated_decay          41256972               human
#> 9           protein_coding          41197695               human
#> 10          protein_coding          41215361               human
#> 11          protein_coding          41215361               human
#> 12          protein_coding          41215377               human
#> 13         retained_intron                NA                <NA>
#> 14          protein_coding          41228505               human
#> 15          protein_coding          41228554               human
#> 16         retained_intron                NA                <NA>
#> 17          non_stop_decay          41245587               human
#> 18          protein_coding          41245601               human
#> 19          protein_coding          41245603               human
#> 20 nonsense_mediated_decay          41262552               human
#> 21          protein_coding          41246129               human
#> 22          protein_coding          41246129               human
#> 23          protein_coding          41246187               human
#> 24          protein_coding          41247863               human
#> 25 nonsense_mediated_decay          41256972               human
#> 26          protein_coding          41256206               human
#> 27          protein_coding          41197695               human
#> 28          protein_coding          41197695               human
#> 29          protein_coding          41197695               human
#> 30          protein_coding          41197695               human
#> 31          protein_coding          41197695               human
#>    Translation.end Translation.object_type Translation.length
#> 1         41276113             Translation               1863
#> 2         41276113             Translation                699
#> 3         41277202             Translation                173
#> 4         41226495             Translation                354
#> 5         41202109             Translation                 96
#> 6         41258543             Translation               1816
#> 7         41276113             Translation               1884
#> 8         41276113             Translation                 63
#> 9         41276113             Translation                759
#> 10        41256933             Translation                498
#> 11        41276113             Translation                623
#> 12        41258543             Translation                572
#> 13              NA                    <NA>                 NA
#> 14        41256933             Translation                266
#> 15        41243841             Translation                242
#> 16              NA                    <NA>                 NA
#> 17        41247883             Translation                437
#> 18        41276113             Translation                649
#> 19        41276113             Translation                622
#> 20        41276113             Translation                 59
#> 21        41246659             Translation                177
#> 22        41276113             Translation                473
#> 23        41256908             Translation                319
#> 24        41276113             Translation                222
#> 25        41276113             Translation                 63
#> 26        41276113             Translation                 98
#> 27        41276113             Translation               1598
#> 28        41276113             Translation                721
#> 29        41276113             Translation               1624
#> 30        41276113             Translation                680
#> 31        41246659             Translation               1567
#>    Translation.version  Translation.id Translation.Parent Translation.db_type
#> 1                    3 ENSP00000350283    ENST00000357654                core
#> 2                    1 ENSP00000417148    ENST00000468300                core
#> 3                    1 ENSP00000465818    ENST00000586385                core
#> 4                    1 ENSP00000467329    ENST00000591534                core
#> 5                    1 ENSP00000465347    ENST00000591849                core
#> 6                    1 ENSP00000418775    ENST00000493795                core
#> 7                    2 ENSP00000418960    ENST00000471181                core
#> 8                    1 ENSP00000418548    ENST00000461221                core
#> 9                    2 ENSP00000420705    ENST00000491747                core
#> 10                   1 ENSP00000419481    ENST00000484087                core
#> 11                   1 ENSP00000420412    ENST00000478531                core
#> 12                   1 ENSP00000418819    ENST00000493919                core
#> 13                  NA            <NA>               <NA>                <NA>
#> 14                   1 ENSP00000418212    ENST00000487825                core
#> 15                   1 ENSP00000417241    ENST00000461574                core
#> 16                  NA            <NA>               <NA>                <NA>
#> 17                   3 ENSP00000397145    ENST00000412061                core
#> 18                   1 ENSP00000419274    ENST00000470026                core
#> 19                   1 ENSP00000419988    ENST00000477152                core
#> 20                   1 ENSP00000420253    ENST00000492859                core
#> 21                   1 ENSP00000418986    ENST00000497488                core
#> 22                   1 ENSP00000419103    ENST00000494123                core
#> 23                   1 ENSP00000420201    ENST00000473961                core
#> 24                   1 ENSP00000417554    ENST00000476777                core
#> 25                   1 ENSP00000417988    ENST00000461798                core
#> 26                   1 ENSP00000420781    ENST00000489037                core
#> 27                   6 ENSP00000326002    ENST00000354071                core
#> 28                   5 ENSP00000312236    ENST00000352993                core
#> 29                   4 ENSP00000246907    ENST00000346315                core
#> 30                   3 ENSP00000338007    ENST00000351666                core
#> 31                   4 ENSP00000310938    ENST00000309486                core
#>            source length    start
#> 1  ensembl_havana   7094 41196312
#> 2  ensembl_havana   3273 41196822
#> 3          havana    781 41197580
#> 4          havana   1282 41197580
#> 5          havana    563 41197580
#> 6  ensembl_havana   5732 41197646
#> 7  ensembl_havana   5936 41197646
#> 8          havana   5693 41197695
#> 9          havana   2379 41197695
#> 10         havana   1495 41215361
#> 11         havana   1972 41215361
#> 12         havana   1948 41215377
#> 13         havana    561 41219291
#> 14         havana    800 41228505
#> 15         havana    726 41228554
#> 16         havana   4497 41243115
#> 17         havana   1312 41245587
#> 18         havana   2108 41245601
#> 19         havana   1980 41245603
#> 20         havana   1584 41246129
#> 21         havana    779 41246129
#> 22         havana   1612 41246129
#> 23         havana    958 41246187
#> 24         havana    769 41247863
#> 25         havana    582 41251848
#> 26         havana    455 41256206
#> 27        ensembl   6411 41196313
#> 28        ensembl   3780 41196313
#> 29        ensembl   6451 41196313
#> 30        ensembl   3444 41196313
#> 31        ensembl   7114 41196313
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
#> 
#> $strand
#> [1] -1
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $version
#> [1] 15
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $db_type
#> [1] "core"
#> 
#> $start
#> [1] 41196312
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
```
