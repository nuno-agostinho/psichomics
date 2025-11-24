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
#> $display_name
#> [1] "BRCA1"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $start
#> [1] 41196312
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
#> $species
#> [1] "human"
#> 
#> $Transcript
#>             Parent         Exon is_canonical display_name object_type
#> 1  ENSG00000012048 c(412773....            0    BRCA1-001  Transcript
#> 2  ENSG00000012048 c("17", ....            0    BRCA1-007  Transcript
#> 3  ENSG00000012048 c("17", ....            0    BRCA1-023  Transcript
#> 4  ENSG00000012048 c("human....            0    BRCA1-024  Transcript
#> 5  ENSG00000012048 c("core"....            0    BRCA1-025  Transcript
#> 6  ENSG00000012048 c("GRCh3....            0    BRCA1-006  Transcript
#> 7  ENSG00000012048 c(1, 1, ....            1    BRCA1-005  Transcript
#> 8  ENSG00000012048 c(412773....            0    BRCA1-010  Transcript
#> 9  ENSG00000012048 c(412773....            0    BRCA1-014  Transcript
#> 10 ENSG00000012048 c(412568....            0    BRCA1-015  Transcript
#> 11 ENSG00000012048 c("GRCh3....            0    BRCA1-009  Transcript
#> 12 ENSG00000012048 c(412774....            0    BRCA1-008  Transcript
#> 13 ENSG00000012048 c(412230....            0    BRCA1-021  Transcript
#> 14 ENSG00000012048 c(412568....            0    BRCA1-019  Transcript
#> 15 ENSG00000012048 c("GRCh3....            0    BRCA1-022  Transcript
#> 16 ENSG00000012048 c("core"....            0    BRCA1-012  Transcript
#> 17 ENSG00000012048 c(-1, -1....            0    BRCA1-026  Transcript
#> 18 ENSG00000012048 c(-1, -1....            0    BRCA1-011  Transcript
#> 19 ENSG00000012048 c(1, 1, ....            0    BRCA1-004  Transcript
#> 20 ENSG00000012048 c(1, 1, ....            0    BRCA1-002  Transcript
#> 21 ENSG00000012048 c(-1, -1....            0    BRCA1-003  Transcript
#> 22 ENSG00000012048 c(412774....            0    BRCA1-013  Transcript
#> 23 ENSG00000012048 c(412569....            0    BRCA1-018  Transcript
#> 24 ENSG00000012048 c("GRCh3....            0    BRCA1-017  Transcript
#> 25 ENSG00000012048 c("GRCh3....            0    BRCA1-020  Transcript
#> 26 ENSG00000012048 c(412771....            0    BRCA1-016  Transcript
#> 27 ENSG00000012048 c("core"....            0    BRCA1-205  Transcript
#> 28 ENSG00000012048 c("GRCh3....            0    BRCA1-204  Transcript
#> 29 ENSG00000012048 c(412774....            0    BRCA1-202  Transcript
#> 30 ENSG00000012048 c("human....            0    BRCA1-203  Transcript
#> 31 ENSG00000012048 c(-1, -1....            0    BRCA1-201  Transcript
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
#>    Translation.db_type  Translation.id Translation.object_type
#> 1                 core ENSP00000350283             Translation
#> 2                 core ENSP00000417148             Translation
#> 3                 core ENSP00000465818             Translation
#> 4                 core ENSP00000467329             Translation
#> 5                 core ENSP00000465347             Translation
#> 6                 core ENSP00000418775             Translation
#> 7                 core ENSP00000418960             Translation
#> 8                 core ENSP00000418548             Translation
#> 9                 core ENSP00000420705             Translation
#> 10                core ENSP00000419481             Translation
#> 11                core ENSP00000420412             Translation
#> 12                core ENSP00000418819             Translation
#> 13                <NA>            <NA>                    <NA>
#> 14                core ENSP00000418212             Translation
#> 15                core ENSP00000417241             Translation
#> 16                <NA>            <NA>                    <NA>
#> 17                core ENSP00000397145             Translation
#> 18                core ENSP00000419274             Translation
#> 19                core ENSP00000419988             Translation
#> 20                core ENSP00000420253             Translation
#> 21                core ENSP00000418986             Translation
#> 22                core ENSP00000419103             Translation
#> 23                core ENSP00000420201             Translation
#> 24                core ENSP00000417554             Translation
#> 25                core ENSP00000417988             Translation
#> 26                core ENSP00000420781             Translation
#> 27                core ENSP00000326002             Translation
#> 28                core ENSP00000312236             Translation
#> 29                core ENSP00000246907             Translation
#> 30                core ENSP00000338007             Translation
#> 31                core ENSP00000310938             Translation
#>    Translation.species Translation.Parent    start                 biotype
#> 1                human    ENST00000357654 41196312          protein_coding
#> 2                human    ENST00000468300 41196822          protein_coding
#> 3                human    ENST00000586385 41197580          protein_coding
#> 4                human    ENST00000591534 41197580          protein_coding
#> 5                human    ENST00000591849 41197580          protein_coding
#> 6                human    ENST00000493795 41197646          protein_coding
#> 7                human    ENST00000471181 41197646          protein_coding
#> 8                human    ENST00000461221 41197695 nonsense_mediated_decay
#> 9                human    ENST00000491747 41197695          protein_coding
#> 10               human    ENST00000484087 41215361          protein_coding
#> 11               human    ENST00000478531 41215361          protein_coding
#> 12               human    ENST00000493919 41215377          protein_coding
#> 13                <NA>               <NA> 41219291         retained_intron
#> 14               human    ENST00000487825 41228505          protein_coding
#> 15               human    ENST00000461574 41228554          protein_coding
#> 16                <NA>               <NA> 41243115         retained_intron
#> 17               human    ENST00000412061 41245587          non_stop_decay
#> 18               human    ENST00000470026 41245601          protein_coding
#> 19               human    ENST00000477152 41245603          protein_coding
#> 20               human    ENST00000492859 41246129 nonsense_mediated_decay
#> 21               human    ENST00000497488 41246129          protein_coding
#> 22               human    ENST00000494123 41246129          protein_coding
#> 23               human    ENST00000473961 41246187          protein_coding
#> 24               human    ENST00000476777 41247863          protein_coding
#> 25               human    ENST00000461798 41251848 nonsense_mediated_decay
#> 26               human    ENST00000489037 41256206          protein_coding
#> 27               human    ENST00000354071 41196313          protein_coding
#> 28               human    ENST00000352993 41196313          protein_coding
#> 29               human    ENST00000346315 41196313          protein_coding
#> 30               human    ENST00000351666 41196313          protein_coding
#> 31               human    ENST00000309486 41196313          protein_coding
#>    length version species seq_region_name assembly_name strand gencode_primary
#> 1    7094       3   human              17        GRCh37     -1               0
#> 2    3273       1   human              17        GRCh37     -1               0
#> 3     781       1   human              17        GRCh37     -1               0
#> 4    1282       1   human              17        GRCh37     -1               0
#> 5     563       1   human              17        GRCh37     -1               0
#> 6    5732       1   human              17        GRCh37     -1               0
#> 7    5936       2   human              17        GRCh37     -1               0
#> 8    5693       1   human              17        GRCh37     -1               0
#> 9    2379       2   human              17        GRCh37     -1               0
#> 10   1495       1   human              17        GRCh37     -1               0
#> 11   1972       1   human              17        GRCh37     -1               0
#> 12   1948       1   human              17        GRCh37     -1               0
#> 13    561       1   human              17        GRCh37     -1               0
#> 14    800       1   human              17        GRCh37     -1               0
#> 15    726       1   human              17        GRCh37     -1               0
#> 16   4497       1   human              17        GRCh37     -1               0
#> 17   1312       3   human              17        GRCh37     -1               0
#> 18   2108       1   human              17        GRCh37     -1               0
#> 19   1980       1   human              17        GRCh37     -1               0
#> 20   1584       1   human              17        GRCh37     -1               0
#> 21    779       1   human              17        GRCh37     -1               0
#> 22   1612       1   human              17        GRCh37     -1               0
#> 23    958       1   human              17        GRCh37     -1               0
#> 24    769       1   human              17        GRCh37     -1               0
#> 25    582       1   human              17        GRCh37     -1               0
#> 26    455       1   human              17        GRCh37     -1               0
#> 27   6411       3   human              17        GRCh37     -1               0
#> 28   3780       3   human              17        GRCh37     -1               0
#> 29   6451       3   human              17        GRCh37     -1               0
#> 30   3444       3   human              17        GRCh37     -1               0
#> 31   7114       4   human              17        GRCh37     -1               0
#>                 id db_type                logic_name      end         source
#> 1  ENST00000357654    core ensembl_havana_transcript 41277387 ensembl_havana
#> 2  ENST00000468300    core ensembl_havana_transcript 41277468 ensembl_havana
#> 3  ENST00000586385    core    havana_homo_sapiens_37 41277346         havana
#> 4  ENST00000591534    core    havana_homo_sapiens_37 41277346         havana
#> 5  ENST00000591849    core    havana_homo_sapiens_37 41277346         havana
#> 6  ENST00000493795    core ensembl_havana_transcript 41277419 ensembl_havana
#> 7  ENST00000471181    core ensembl_havana_transcript 41277500 ensembl_havana
#> 8  ENST00000461221    core    havana_homo_sapiens_37 41277305         havana
#> 9  ENST00000491747    core    havana_homo_sapiens_37 41277373         havana
#> 10 ENST00000484087    core    havana_homo_sapiens_37 41256933         havana
#> 11 ENST00000478531    core    havana_homo_sapiens_37 41277376         havana
#> 12 ENST00000493919    core    havana_homo_sapiens_37 41277419         havana
#> 13 ENST00000472490    core    havana_homo_sapiens_37 41223083         havana
#> 14 ENST00000487825    core    havana_homo_sapiens_37 41256933         havana
#> 15 ENST00000461574    core    havana_homo_sapiens_37 41243841         havana
#> 16 ENST00000467274    core    havana_homo_sapiens_37 41277332         havana
#> 17 ENST00000412061    core    havana_homo_sapiens_37 41247883         havana
#> 18 ENST00000470026    core    havana_homo_sapiens_37 41277340         havana
#> 19 ENST00000477152    core    havana_homo_sapiens_37 41277381         havana
#> 20 ENST00000492859    core    havana_homo_sapiens_37 41277317         havana
#> 21 ENST00000497488    core    havana_homo_sapiens_37 41277317         havana
#> 22 ENST00000494123    core    havana_homo_sapiens_37 41277467         havana
#> 23 ENST00000473961    core    havana_homo_sapiens_37 41256908         havana
#> 24 ENST00000476777    core    havana_homo_sapiens_37 41277370         havana
#> 25 ENST00000461798    core    havana_homo_sapiens_37 41277387         havana
#> 26 ENST00000489037    core    havana_homo_sapiens_37 41277338         havana
#> 27 ENST00000354071    core   ensembl_homo_sapiens_37 41277500        ensembl
#> 28 ENST00000352993    core   ensembl_homo_sapiens_37 41277500        ensembl
#> 29 ENST00000346315    core   ensembl_homo_sapiens_37 41277468        ensembl
#> 30 ENST00000351666    core   ensembl_homo_sapiens_37 41276132        ensembl
#> 31 ENST00000309486    core   ensembl_homo_sapiens_37 41277468        ensembl
#> 
#> $seq_region_name
#> [1] "17"
#> 
#> $db_type
#> [1] "core"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $strand
#> [1] -1
#> 
#> $end
#> [1] 41277500
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
queryEnsemblByGene("ENSG00000139618")
#> $start
#> [1] 32889611
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $version
#> [1] 10
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
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $db_type
#> [1] "core"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $species
#> [1] "homo_sapiens"
#> 
#> $Transcript
#>           source db_type object_type          Parent      end length
#> 1 ensembl_havana    core  Transcript ENSG00000139618 32973347  10930
#> 2         havana    core  Transcript ENSG00000139618 32907428   2011
#> 3         havana    core  Transcript ENSG00000139618 32953632    495
#> 4         havana    core  Transcript ENSG00000139618 32972409    842
#> 5         havana    core  Transcript ENSG00000139618 32972585    523
#> 6        ensembl    core  Transcript ENSG00000139618 32973805  10984
#>   display_name version    start assembly_name gencode_primary
#> 1    BRCA2-001       3 32889611        GRCh37               0
#> 2    BRCA2-003       2 32889642        GRCh37               0
#> 3    BRCA2-005       1 32945108        GRCh37               0
#> 4    BRCA2-002       1 32953977        GRCh37               0
#> 5    BRCA2-006       1 32970946        GRCh37               0
#> 6    BRCA2-201       1 32889617        GRCh37               0
#>                   biotype                logic_name              id
#> 1          protein_coding ensembl_havana_transcript ENST00000380152
#> 2          protein_coding    havana_homo_sapiens_37 ENST00000530893
#> 3 nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000528762
#> 4 nonsense_mediated_decay    havana_homo_sapiens_37 ENST00000470094
#> 5         retained_intron    havana_homo_sapiens_37 ENST00000533776
#> 6          protein_coding   ensembl_homo_sapiens_37 ENST00000544455
#>   Translation.length Translation.end Translation.version Translation.species
#> 1               3418        32972907                   3        homo_sapiens
#> 2                481        32907428                   2        homo_sapiens
#> 3                 64        32950807                   1        homo_sapiens
#> 4                186        32970229                   1        homo_sapiens
#> 5                 NA              NA                  NA                <NA>
#> 6               3418        32972907                   1        homo_sapiens
#>   Translation.start Translation.db_type Translation.object_type
#> 1          32890598                core             Translation
#> 2          32899266                core             Translation
#> 3          32945108                core             Translation
#> 4          32953977                core             Translation
#> 5                NA                <NA>                    <NA>
#> 6          32890598                core             Translation
#>   Translation.Parent  Translation.id strand         Exon seq_region_name
#> 1    ENST00000380152 ENSP00000369497      1 c(1, 1, ....              13
#> 2    ENST00000530893 ENSP00000435699      1 c("ENSE0....              13
#> 3    ENST00000528762 ENSP00000433168      1 c("13", ....              13
#> 4    ENST00000470094 ENSP00000434898      1 c("Exon"....              13
#> 5               <NA>            <NA>      1 c(329709....              13
#> 6    ENST00000544455 ENSP00000439902      1 c(328898....              13
#>   is_canonical      species
#> 1            0 homo_sapiens
#> 2            0 homo_sapiens
#> 3            0 homo_sapiens
#> 4            0 homo_sapiens
#> 5            0 homo_sapiens
#> 6            1 homo_sapiens
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $strand
#> [1] 1
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $Transcript
#>       start Translation.Parent  Translation.id Translation.end
#> 1  41196312    ENST00000357654 ENSP00000350283        41276113
#> 2  41196822    ENST00000468300 ENSP00000417148        41276113
#> 3  41197580    ENST00000586385 ENSP00000465818        41277202
#> 4  41197580    ENST00000591534 ENSP00000467329        41226495
#> 5  41197580    ENST00000591849 ENSP00000465347        41202109
#> 6  41197646    ENST00000493795 ENSP00000418775        41258543
#> 7  41197646    ENST00000471181 ENSP00000418960        41276113
#> 8  41197695    ENST00000461221 ENSP00000418548        41276113
#> 9  41197695    ENST00000491747 ENSP00000420705        41276113
#> 10 41215361    ENST00000484087 ENSP00000419481        41256933
#> 11 41215361    ENST00000478531 ENSP00000420412        41276113
#> 12 41215377    ENST00000493919 ENSP00000418819        41258543
#> 13 41219291               <NA>            <NA>              NA
#> 14 41228505    ENST00000487825 ENSP00000418212        41256933
#> 15 41228554    ENST00000461574 ENSP00000417241        41243841
#> 16 41243115               <NA>            <NA>              NA
#> 17 41245587    ENST00000412061 ENSP00000397145        41247883
#> 18 41245601    ENST00000470026 ENSP00000419274        41276113
#> 19 41245603    ENST00000477152 ENSP00000419988        41276113
#> 20 41246129    ENST00000492859 ENSP00000420253        41276113
#> 21 41246129    ENST00000497488 ENSP00000418986        41246659
#> 22 41246129    ENST00000494123 ENSP00000419103        41276113
#> 23 41246187    ENST00000473961 ENSP00000420201        41256908
#> 24 41247863    ENST00000476777 ENSP00000417554        41276113
#> 25 41251848    ENST00000461798 ENSP00000417988        41276113
#> 26 41256206    ENST00000489037 ENSP00000420781        41276113
#> 27 41196313    ENST00000354071 ENSP00000326002        41276113
#> 28 41196313    ENST00000352993 ENSP00000312236        41276113
#> 29 41196313    ENST00000346315 ENSP00000246907        41276113
#> 30 41196313    ENST00000351666 ENSP00000338007        41276113
#> 31 41196313    ENST00000309486 ENSP00000310938        41246659
#>    Translation.version Translation.object_type Translation.db_type
#> 1                    3             Translation                core
#> 2                    1             Translation                core
#> 3                    1             Translation                core
#> 4                    1             Translation                core
#> 5                    1             Translation                core
#> 6                    1             Translation                core
#> 7                    2             Translation                core
#> 8                    1             Translation                core
#> 9                    2             Translation                core
#> 10                   1             Translation                core
#> 11                   1             Translation                core
#> 12                   1             Translation                core
#> 13                  NA                    <NA>                <NA>
#> 14                   1             Translation                core
#> 15                   1             Translation                core
#> 16                  NA                    <NA>                <NA>
#> 17                   3             Translation                core
#> 18                   1             Translation                core
#> 19                   1             Translation                core
#> 20                   1             Translation                core
#> 21                   1             Translation                core
#> 22                   1             Translation                core
#> 23                   1             Translation                core
#> 24                   1             Translation                core
#> 25                   1             Translation                core
#> 26                   1             Translation                core
#> 27                   6             Translation                core
#> 28                   5             Translation                core
#> 29                   4             Translation                core
#> 30                   3             Translation                core
#> 31                   4             Translation                core
#>    Translation.length Translation.start Translation.species display_name
#> 1                1863          41197695               human    BRCA1-001
#> 2                 699          41197801               human    BRCA1-007
#> 3                 173          41197695               human    BRCA1-023
#> 4                 354          41197695               human    BRCA1-024
#> 5                  96          41197695               human    BRCA1-025
#> 6                1816          41197695               human    BRCA1-006
#> 7                1884          41197695               human    BRCA1-005
#> 8                  63          41256972               human    BRCA1-010
#> 9                 759          41197695               human    BRCA1-014
#> 10                498          41215361               human    BRCA1-015
#> 11                623          41215361               human    BRCA1-009
#> 12                572          41215377               human    BRCA1-008
#> 13                 NA                NA                <NA>    BRCA1-021
#> 14                266          41228505               human    BRCA1-019
#> 15                242          41228554               human    BRCA1-022
#> 16                 NA                NA                <NA>    BRCA1-012
#> 17                437          41245587               human    BRCA1-026
#> 18                649          41245601               human    BRCA1-011
#> 19                622          41245603               human    BRCA1-004
#> 20                 59          41262552               human    BRCA1-002
#> 21                177          41246129               human    BRCA1-003
#> 22                473          41246129               human    BRCA1-013
#> 23                319          41246187               human    BRCA1-018
#> 24                222          41247863               human    BRCA1-017
#> 25                 63          41256972               human    BRCA1-020
#> 26                 98          41256206               human    BRCA1-016
#> 27               1598          41197695               human    BRCA1-205
#> 28                721          41197695               human    BRCA1-204
#> 29               1624          41197695               human    BRCA1-202
#> 30                680          41197695               human    BRCA1-203
#> 31               1567          41197695               human    BRCA1-201
#>            source length assembly_name version                logic_name
#> 1  ensembl_havana   7094        GRCh37       3 ensembl_havana_transcript
#> 2  ensembl_havana   3273        GRCh37       1 ensembl_havana_transcript
#> 3          havana    781        GRCh37       1    havana_homo_sapiens_37
#> 4          havana   1282        GRCh37       1    havana_homo_sapiens_37
#> 5          havana    563        GRCh37       1    havana_homo_sapiens_37
#> 6  ensembl_havana   5732        GRCh37       1 ensembl_havana_transcript
#> 7  ensembl_havana   5936        GRCh37       2 ensembl_havana_transcript
#> 8          havana   5693        GRCh37       1    havana_homo_sapiens_37
#> 9          havana   2379        GRCh37       2    havana_homo_sapiens_37
#> 10         havana   1495        GRCh37       1    havana_homo_sapiens_37
#> 11         havana   1972        GRCh37       1    havana_homo_sapiens_37
#> 12         havana   1948        GRCh37       1    havana_homo_sapiens_37
#> 13         havana    561        GRCh37       1    havana_homo_sapiens_37
#> 14         havana    800        GRCh37       1    havana_homo_sapiens_37
#> 15         havana    726        GRCh37       1    havana_homo_sapiens_37
#> 16         havana   4497        GRCh37       1    havana_homo_sapiens_37
#> 17         havana   1312        GRCh37       3    havana_homo_sapiens_37
#> 18         havana   2108        GRCh37       1    havana_homo_sapiens_37
#> 19         havana   1980        GRCh37       1    havana_homo_sapiens_37
#> 20         havana   1584        GRCh37       1    havana_homo_sapiens_37
#> 21         havana    779        GRCh37       1    havana_homo_sapiens_37
#> 22         havana   1612        GRCh37       1    havana_homo_sapiens_37
#> 23         havana    958        GRCh37       1    havana_homo_sapiens_37
#> 24         havana    769        GRCh37       1    havana_homo_sapiens_37
#> 25         havana    582        GRCh37       1    havana_homo_sapiens_37
#> 26         havana    455        GRCh37       1    havana_homo_sapiens_37
#> 27        ensembl   6411        GRCh37       3   ensembl_homo_sapiens_37
#> 28        ensembl   3780        GRCh37       3   ensembl_homo_sapiens_37
#> 29        ensembl   6451        GRCh37       3   ensembl_homo_sapiens_37
#> 30        ensembl   3444        GRCh37       3   ensembl_homo_sapiens_37
#> 31        ensembl   7114        GRCh37       4   ensembl_homo_sapiens_37
#>    object_type                 biotype gencode_primary              id
#> 1   Transcript          protein_coding               0 ENST00000357654
#> 2   Transcript          protein_coding               0 ENST00000468300
#> 3   Transcript          protein_coding               0 ENST00000586385
#> 4   Transcript          protein_coding               0 ENST00000591534
#> 5   Transcript          protein_coding               0 ENST00000591849
#> 6   Transcript          protein_coding               0 ENST00000493795
#> 7   Transcript          protein_coding               0 ENST00000471181
#> 8   Transcript nonsense_mediated_decay               0 ENST00000461221
#> 9   Transcript          protein_coding               0 ENST00000491747
#> 10  Transcript          protein_coding               0 ENST00000484087
#> 11  Transcript          protein_coding               0 ENST00000478531
#> 12  Transcript          protein_coding               0 ENST00000493919
#> 13  Transcript         retained_intron               0 ENST00000472490
#> 14  Transcript          protein_coding               0 ENST00000487825
#> 15  Transcript          protein_coding               0 ENST00000461574
#> 16  Transcript         retained_intron               0 ENST00000467274
#> 17  Transcript          non_stop_decay               0 ENST00000412061
#> 18  Transcript          protein_coding               0 ENST00000470026
#> 19  Transcript          protein_coding               0 ENST00000477152
#> 20  Transcript nonsense_mediated_decay               0 ENST00000492859
#> 21  Transcript          protein_coding               0 ENST00000497488
#> 22  Transcript          protein_coding               0 ENST00000494123
#> 23  Transcript          protein_coding               0 ENST00000473961
#> 24  Transcript          protein_coding               0 ENST00000476777
#> 25  Transcript nonsense_mediated_decay               0 ENST00000461798
#> 26  Transcript          protein_coding               0 ENST00000489037
#> 27  Transcript          protein_coding               0 ENST00000354071
#> 28  Transcript          protein_coding               0 ENST00000352993
#> 29  Transcript          protein_coding               0 ENST00000346315
#> 30  Transcript          protein_coding               0 ENST00000351666
#> 31  Transcript          protein_coding               0 ENST00000309486
#>    is_canonical species         Exon strand db_type          Parent
#> 1             0   human c(-1, -1....     -1    core ENSG00000012048
#> 2             0   human c(1, 1, ....     -1    core ENSG00000012048
#> 3             0   human c("Exon"....     -1    core ENSG00000012048
#> 4             0   human c(-1, -1....     -1    core ENSG00000012048
#> 5             0   human c(-1, -1....     -1    core ENSG00000012048
#> 6             0   human c(-1, -1....     -1    core ENSG00000012048
#> 7             1   human c(-1, -1....     -1    core ENSG00000012048
#> 8             0   human c(412773....     -1    core ENSG00000012048
#> 9             0   human c(412772....     -1    core ENSG00000012048
#> 10            0   human c(-1, -1....     -1    core ENSG00000012048
#> 11            0   human c(1, 1, ....     -1    core ENSG00000012048
#> 12            0   human c("human....     -1    core ENSG00000012048
#> 13            0   human c(412229....     -1    core ENSG00000012048
#> 14            0   human c(412568....     -1    core ENSG00000012048
#> 15            0   human c(412434....     -1    core ENSG00000012048
#> 16            0   human c(-1, -1....     -1    core ENSG00000012048
#> 17            0   human c("GRCh3....     -1    core ENSG00000012048
#> 18            0   human c("17", ....     -1    core ENSG00000012048
#> 19            0   human c("core"....     -1    core ENSG00000012048
#> 20            0   human c("Exon"....     -1    core ENSG00000012048
#> 21            0   human c(1, 1),....     -1    core ENSG00000012048
#> 22            0   human c("17", ....     -1    core ENSG00000012048
#> 23            0   human c(412568....     -1    core ENSG00000012048
#> 24            0   human c(-1, -1....     -1    core ENSG00000012048
#> 25            0   human c(412772....     -1    core ENSG00000012048
#> 26            0   human c(412771....     -1    core ENSG00000012048
#> 27            0   human c(-1, -1....     -1    core ENSG00000012048
#> 28            0   human c("core"....     -1    core ENSG00000012048
#> 29            0   human c(1, 1, ....     -1    core ENSG00000012048
#> 30            0   human c("17", ....     -1    core ENSG00000012048
#> 31            0   human c("17", ....     -1    core ENSG00000012048
#>    seq_region_name      end
#> 1               17 41277387
#> 2               17 41277468
#> 3               17 41277346
#> 4               17 41277346
#> 5               17 41277346
#> 6               17 41277419
#> 7               17 41277500
#> 8               17 41277305
#> 9               17 41277373
#> 10              17 41256933
#> 11              17 41277376
#> 12              17 41277419
#> 13              17 41223083
#> 14              17 41256933
#> 15              17 41243841
#> 16              17 41277332
#> 17              17 41247883
#> 18              17 41277340
#> 19              17 41277381
#> 20              17 41277317
#> 21              17 41277317
#> 22              17 41277467
#> 23              17 41256908
#> 24              17 41277370
#> 25              17 41277387
#> 26              17 41277338
#> 27              17 41277500
#> 28              17 41277500
#> 29              17 41277468
#> 30              17 41276132
#> 31              17 41277468
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $version
#> [1] 15
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $object_type
#> [1] "Gene"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $start
#> [1] 41196312
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $seq_region_name
#> [1] "17"
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
#> $strand
#> [1] -1
#> 
#> $species
#> [1] "human"
#> 
```
