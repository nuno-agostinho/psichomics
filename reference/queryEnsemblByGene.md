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
#> $db_type
#> [1] "core"
#> 
#> $object_type
#> [1] "Gene"
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
#> $version
#> [1] 15
#> 
#> $Transcript
#>            source is_canonical          Parent                 biotype species
#> 1  ensembl_havana            0 ENSG00000012048          protein_coding   human
#> 2  ensembl_havana            0 ENSG00000012048          protein_coding   human
#> 3          havana            0 ENSG00000012048          protein_coding   human
#> 4          havana            0 ENSG00000012048          protein_coding   human
#> 5          havana            0 ENSG00000012048          protein_coding   human
#> 6  ensembl_havana            0 ENSG00000012048          protein_coding   human
#> 7  ensembl_havana            1 ENSG00000012048          protein_coding   human
#> 8          havana            0 ENSG00000012048 nonsense_mediated_decay   human
#> 9          havana            0 ENSG00000012048          protein_coding   human
#> 10         havana            0 ENSG00000012048          protein_coding   human
#> 11         havana            0 ENSG00000012048          protein_coding   human
#> 12         havana            0 ENSG00000012048          protein_coding   human
#> 13         havana            0 ENSG00000012048         retained_intron   human
#> 14         havana            0 ENSG00000012048          protein_coding   human
#> 15         havana            0 ENSG00000012048          protein_coding   human
#> 16         havana            0 ENSG00000012048         retained_intron   human
#> 17         havana            0 ENSG00000012048          non_stop_decay   human
#> 18         havana            0 ENSG00000012048          protein_coding   human
#> 19         havana            0 ENSG00000012048          protein_coding   human
#> 20         havana            0 ENSG00000012048 nonsense_mediated_decay   human
#> 21         havana            0 ENSG00000012048          protein_coding   human
#> 22         havana            0 ENSG00000012048          protein_coding   human
#> 23         havana            0 ENSG00000012048          protein_coding   human
#> 24         havana            0 ENSG00000012048          protein_coding   human
#> 25         havana            0 ENSG00000012048 nonsense_mediated_decay   human
#> 26         havana            0 ENSG00000012048          protein_coding   human
#> 27        ensembl            0 ENSG00000012048          protein_coding   human
#> 28        ensembl            0 ENSG00000012048          protein_coding   human
#> 29        ensembl            0 ENSG00000012048          protein_coding   human
#> 30        ensembl            0 ENSG00000012048          protein_coding   human
#> 31        ensembl            0 ENSG00000012048          protein_coding   human
#>    seq_region_name                logic_name strand      end object_type
#> 1               17 ensembl_havana_transcript     -1 41277387  Transcript
#> 2               17 ensembl_havana_transcript     -1 41277468  Transcript
#> 3               17    havana_homo_sapiens_37     -1 41277346  Transcript
#> 4               17    havana_homo_sapiens_37     -1 41277346  Transcript
#> 5               17    havana_homo_sapiens_37     -1 41277346  Transcript
#> 6               17 ensembl_havana_transcript     -1 41277419  Transcript
#> 7               17 ensembl_havana_transcript     -1 41277500  Transcript
#> 8               17    havana_homo_sapiens_37     -1 41277305  Transcript
#> 9               17    havana_homo_sapiens_37     -1 41277373  Transcript
#> 10              17    havana_homo_sapiens_37     -1 41256933  Transcript
#> 11              17    havana_homo_sapiens_37     -1 41277376  Transcript
#> 12              17    havana_homo_sapiens_37     -1 41277419  Transcript
#> 13              17    havana_homo_sapiens_37     -1 41223083  Transcript
#> 14              17    havana_homo_sapiens_37     -1 41256933  Transcript
#> 15              17    havana_homo_sapiens_37     -1 41243841  Transcript
#> 16              17    havana_homo_sapiens_37     -1 41277332  Transcript
#> 17              17    havana_homo_sapiens_37     -1 41247883  Transcript
#> 18              17    havana_homo_sapiens_37     -1 41277340  Transcript
#> 19              17    havana_homo_sapiens_37     -1 41277381  Transcript
#> 20              17    havana_homo_sapiens_37     -1 41277317  Transcript
#> 21              17    havana_homo_sapiens_37     -1 41277317  Transcript
#> 22              17    havana_homo_sapiens_37     -1 41277467  Transcript
#> 23              17    havana_homo_sapiens_37     -1 41256908  Transcript
#> 24              17    havana_homo_sapiens_37     -1 41277370  Transcript
#> 25              17    havana_homo_sapiens_37     -1 41277387  Transcript
#> 26              17    havana_homo_sapiens_37     -1 41277338  Transcript
#> 27              17   ensembl_homo_sapiens_37     -1 41277500  Transcript
#> 28              17   ensembl_homo_sapiens_37     -1 41277500  Transcript
#> 29              17   ensembl_homo_sapiens_37     -1 41277468  Transcript
#> 30              17   ensembl_homo_sapiens_37     -1 41276132  Transcript
#> 31              17   ensembl_homo_sapiens_37     -1 41277468  Transcript
#>            Exon Translation.object_type Translation.db_type Translation.Parent
#> 1  c("core"....             Translation                core    ENST00000357654
#> 2  c("core"....             Translation                core    ENST00000468300
#> 3  c("core"....             Translation                core    ENST00000586385
#> 4  c(412772....             Translation                core    ENST00000591534
#> 5  c(1, 1, ....             Translation                core    ENST00000591849
#> 6  c(1, 1, ....             Translation                core    ENST00000493795
#> 7  c("Exon"....             Translation                core    ENST00000471181
#> 8  c("Exon"....             Translation                core    ENST00000461221
#> 9  c("Exon"....             Translation                core    ENST00000491747
#> 10 c("core"....             Translation                core    ENST00000484087
#> 11 c("core"....             Translation                core    ENST00000478531
#> 12 c(-1, -1....             Translation                core    ENST00000493919
#> 13 c("core"....                    <NA>                <NA>               <NA>
#> 14 c(-1, -1....             Translation                core    ENST00000487825
#> 15 c(-1, -1....             Translation                core    ENST00000461574
#> 16 c(412773....                    <NA>                <NA>               <NA>
#> 17 c(412478....             Translation                core    ENST00000412061
#> 18 c("GRCh3....             Translation                core    ENST00000470026
#> 19 c("core"....             Translation                core    ENST00000477152
#> 20 c("Exon"....             Translation                core    ENST00000492859
#> 21 c("Exon"....             Translation                core    ENST00000497488
#> 22 c(412774....             Translation                core    ENST00000494123
#> 23 c(412569....             Translation                core    ENST00000473961
#> 24 c("core"....             Translation                core    ENST00000476777
#> 25 c(412772....             Translation                core    ENST00000461798
#> 26 c("ENSE0....             Translation                core    ENST00000489037
#> 27 c("core"....             Translation                core    ENST00000354071
#> 28 c(412775....             Translation                core    ENST00000352993
#> 29 c("Exon"....             Translation                core    ENST00000346315
#> 30 c("Exon"....             Translation                core    ENST00000351666
#> 31 c("GRCh3....             Translation                core    ENST00000309486
#>    Translation.length Translation.start Translation.version Translation.species
#> 1                1863          41197695                   3               human
#> 2                 699          41197801                   1               human
#> 3                 173          41197695                   1               human
#> 4                 354          41197695                   1               human
#> 5                  96          41197695                   1               human
#> 6                1816          41197695                   1               human
#> 7                1884          41197695                   2               human
#> 8                  63          41256972                   1               human
#> 9                 759          41197695                   2               human
#> 10                498          41215361                   1               human
#> 11                623          41215361                   1               human
#> 12                572          41215377                   1               human
#> 13                 NA                NA                  NA                <NA>
#> 14                266          41228505                   1               human
#> 15                242          41228554                   1               human
#> 16                 NA                NA                  NA                <NA>
#> 17                437          41245587                   3               human
#> 18                649          41245601                   1               human
#> 19                622          41245603                   1               human
#> 20                 59          41262552                   1               human
#> 21                177          41246129                   1               human
#> 22                473          41246129                   1               human
#> 23                319          41246187                   1               human
#> 24                222          41247863                   1               human
#> 25                 63          41256972                   1               human
#> 26                 98          41256206                   1               human
#> 27               1598          41197695                   6               human
#> 28                721          41197695                   5               human
#> 29               1624          41197695                   4               human
#> 30                680          41197695                   3               human
#> 31               1567          41197695                   4               human
#>     Translation.id Translation.end db_type    start length version
#> 1  ENSP00000350283        41276113    core 41196312   7094       3
#> 2  ENSP00000417148        41276113    core 41196822   3273       1
#> 3  ENSP00000465818        41277202    core 41197580    781       1
#> 4  ENSP00000467329        41226495    core 41197580   1282       1
#> 5  ENSP00000465347        41202109    core 41197580    563       1
#> 6  ENSP00000418775        41258543    core 41197646   5732       1
#> 7  ENSP00000418960        41276113    core 41197646   5936       2
#> 8  ENSP00000418548        41276113    core 41197695   5693       1
#> 9  ENSP00000420705        41276113    core 41197695   2379       2
#> 10 ENSP00000419481        41256933    core 41215361   1495       1
#> 11 ENSP00000420412        41276113    core 41215361   1972       1
#> 12 ENSP00000418819        41258543    core 41215377   1948       1
#> 13            <NA>              NA    core 41219291    561       1
#> 14 ENSP00000418212        41256933    core 41228505    800       1
#> 15 ENSP00000417241        41243841    core 41228554    726       1
#> 16            <NA>              NA    core 41243115   4497       1
#> 17 ENSP00000397145        41247883    core 41245587   1312       3
#> 18 ENSP00000419274        41276113    core 41245601   2108       1
#> 19 ENSP00000419988        41276113    core 41245603   1980       1
#> 20 ENSP00000420253        41276113    core 41246129   1584       1
#> 21 ENSP00000418986        41246659    core 41246129    779       1
#> 22 ENSP00000419103        41276113    core 41246129   1612       1
#> 23 ENSP00000420201        41256908    core 41246187    958       1
#> 24 ENSP00000417554        41276113    core 41247863    769       1
#> 25 ENSP00000417988        41276113    core 41251848    582       1
#> 26 ENSP00000420781        41276113    core 41256206    455       1
#> 27 ENSP00000326002        41276113    core 41196313   6411       3
#> 28 ENSP00000312236        41276113    core 41196313   3780       3
#> 29 ENSP00000246907        41276113    core 41196313   6451       3
#> 30 ENSP00000338007        41276113    core 41196313   3444       3
#> 31 ENSP00000310938        41246659    core 41196313   7114       4
#>    gencode_primary              id display_name assembly_name
#> 1                0 ENST00000357654    BRCA1-001        GRCh37
#> 2                0 ENST00000468300    BRCA1-007        GRCh37
#> 3                0 ENST00000586385    BRCA1-023        GRCh37
#> 4                0 ENST00000591534    BRCA1-024        GRCh37
#> 5                0 ENST00000591849    BRCA1-025        GRCh37
#> 6                0 ENST00000493795    BRCA1-006        GRCh37
#> 7                0 ENST00000471181    BRCA1-005        GRCh37
#> 8                0 ENST00000461221    BRCA1-010        GRCh37
#> 9                0 ENST00000491747    BRCA1-014        GRCh37
#> 10               0 ENST00000484087    BRCA1-015        GRCh37
#> 11               0 ENST00000478531    BRCA1-009        GRCh37
#> 12               0 ENST00000493919    BRCA1-008        GRCh37
#> 13               0 ENST00000472490    BRCA1-021        GRCh37
#> 14               0 ENST00000487825    BRCA1-019        GRCh37
#> 15               0 ENST00000461574    BRCA1-022        GRCh37
#> 16               0 ENST00000467274    BRCA1-012        GRCh37
#> 17               0 ENST00000412061    BRCA1-026        GRCh37
#> 18               0 ENST00000470026    BRCA1-011        GRCh37
#> 19               0 ENST00000477152    BRCA1-004        GRCh37
#> 20               0 ENST00000492859    BRCA1-002        GRCh37
#> 21               0 ENST00000497488    BRCA1-003        GRCh37
#> 22               0 ENST00000494123    BRCA1-013        GRCh37
#> 23               0 ENST00000473961    BRCA1-018        GRCh37
#> 24               0 ENST00000476777    BRCA1-017        GRCh37
#> 25               0 ENST00000461798    BRCA1-020        GRCh37
#> 26               0 ENST00000489037    BRCA1-016        GRCh37
#> 27               0 ENST00000354071    BRCA1-205        GRCh37
#> 28               0 ENST00000352993    BRCA1-204        GRCh37
#> 29               0 ENST00000346315    BRCA1-202        GRCh37
#> 30               0 ENST00000351666    BRCA1-203        GRCh37
#> 31               0 ENST00000309486    BRCA1-201        GRCh37
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
queryEnsemblByGene("ENSG00000139618")
#> $species
#> [1] "homo_sapiens"
#> 
#> $strand
#> [1] 1
#> 
#> $version
#> [1] 10
#> 
#> $canonical_transcript
#> [1] "ENST00000544455.1"
#> 
#> $source
#> [1] "ensembl_havana"
#> 
#> $id
#> [1] "ENSG00000139618"
#> 
#> $seq_region_name
#> [1] "13"
#> 
#> $display_name
#> [1] "BRCA2"
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
#> $description
#> [1] "breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $db_type
#> [1] "core"
#> 
#> $end
#> [1] 32973805
#> 
#> $Transcript
#>   is_canonical      species strand version display_name assembly_name
#> 1            0 homo_sapiens      1       3    BRCA2-001        GRCh37
#> 2            0 homo_sapiens      1       2    BRCA2-003        GRCh37
#> 3            0 homo_sapiens      1       1    BRCA2-005        GRCh37
#> 4            0 homo_sapiens      1       1    BRCA2-002        GRCh37
#> 5            0 homo_sapiens      1       1    BRCA2-006        GRCh37
#> 6            1 homo_sapiens      1       1    BRCA2-201        GRCh37
#>           source          Parent seq_region_name Translation.object_type
#> 1 ensembl_havana ENSG00000139618              13             Translation
#> 2         havana ENSG00000139618              13             Translation
#> 3         havana ENSG00000139618              13             Translation
#> 4         havana ENSG00000139618              13             Translation
#> 5         havana ENSG00000139618              13                    <NA>
#> 6        ensembl ENSG00000139618              13             Translation
#>   Translation.Parent Translation.length  Translation.id Translation.db_type
#> 1    ENST00000380152               3418 ENSP00000369497                core
#> 2    ENST00000530893                481 ENSP00000435699                core
#> 3    ENST00000528762                 64 ENSP00000433168                core
#> 4    ENST00000470094                186 ENSP00000434898                core
#> 5               <NA>                 NA            <NA>                <NA>
#> 6    ENST00000544455               3418 ENSP00000439902                core
#>   Translation.end Translation.start Translation.species Translation.version
#> 1        32972907          32890598        homo_sapiens                   3
#> 2        32907428          32899266        homo_sapiens                   2
#> 3        32950807          32945108        homo_sapiens                   1
#> 4        32970229          32953977        homo_sapiens                   1
#> 5              NA                NA                <NA>                  NA
#> 6        32972907          32890598        homo_sapiens                   1
#>                id                logic_name gencode_primary         Exon
#> 1 ENST00000380152 ensembl_havana_transcript               0 c(4, 1, ....
#> 2 ENST00000530893    havana_homo_sapiens_37               0 c(1, 1, ....
#> 3 ENST00000528762    havana_homo_sapiens_37               0 c(1, 1, ....
#> 4 ENST00000470094    havana_homo_sapiens_37               0 c(1, 1, ....
#> 5 ENST00000533776    havana_homo_sapiens_37               0 c(1, 1),....
#> 6 ENST00000544455   ensembl_homo_sapiens_37               0 c("13", ....
#>        end db_type    start object_type                 biotype length
#> 1 32973347    core 32889611  Transcript          protein_coding  10930
#> 2 32907428    core 32889642  Transcript          protein_coding   2011
#> 3 32953632    core 32945108  Transcript nonsense_mediated_decay    495
#> 4 32972409    core 32953977  Transcript nonsense_mediated_decay    842
#> 5 32972585    core 32970946  Transcript         retained_intron    523
#> 6 32973805    core 32889617  Transcript          protein_coding  10984
#> 
#> $start
#> [1] 32889611
#> 
event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
queryEnsemblByEvent(event, species="human", assembly="hg19")
#> $species
#> [1] "human"
#> 
#> $Transcript
#>    gencode_primary                logic_name                 biotype
#> 1                0 ensembl_havana_transcript          protein_coding
#> 2                0 ensembl_havana_transcript          protein_coding
#> 3                0    havana_homo_sapiens_37          protein_coding
#> 4                0    havana_homo_sapiens_37          protein_coding
#> 5                0    havana_homo_sapiens_37          protein_coding
#> 6                0 ensembl_havana_transcript          protein_coding
#> 7                0 ensembl_havana_transcript          protein_coding
#> 8                0    havana_homo_sapiens_37 nonsense_mediated_decay
#> 9                0    havana_homo_sapiens_37          protein_coding
#> 10               0    havana_homo_sapiens_37          protein_coding
#> 11               0    havana_homo_sapiens_37          protein_coding
#> 12               0    havana_homo_sapiens_37          protein_coding
#> 13               0    havana_homo_sapiens_37         retained_intron
#> 14               0    havana_homo_sapiens_37          protein_coding
#> 15               0    havana_homo_sapiens_37          protein_coding
#> 16               0    havana_homo_sapiens_37         retained_intron
#> 17               0    havana_homo_sapiens_37          non_stop_decay
#> 18               0    havana_homo_sapiens_37          protein_coding
#> 19               0    havana_homo_sapiens_37          protein_coding
#> 20               0    havana_homo_sapiens_37 nonsense_mediated_decay
#> 21               0    havana_homo_sapiens_37          protein_coding
#> 22               0    havana_homo_sapiens_37          protein_coding
#> 23               0    havana_homo_sapiens_37          protein_coding
#> 24               0    havana_homo_sapiens_37          protein_coding
#> 25               0    havana_homo_sapiens_37 nonsense_mediated_decay
#> 26               0    havana_homo_sapiens_37          protein_coding
#> 27               0   ensembl_homo_sapiens_37          protein_coding
#> 28               0   ensembl_homo_sapiens_37          protein_coding
#> 29               0   ensembl_homo_sapiens_37          protein_coding
#> 30               0   ensembl_homo_sapiens_37          protein_coding
#> 31               0   ensembl_homo_sapiens_37          protein_coding
#>                 id Translation.version Translation.end Translation.length
#> 1  ENST00000357654                   3        41276113               1863
#> 2  ENST00000468300                   1        41276113                699
#> 3  ENST00000586385                   1        41277202                173
#> 4  ENST00000591534                   1        41226495                354
#> 5  ENST00000591849                   1        41202109                 96
#> 6  ENST00000493795                   1        41258543               1816
#> 7  ENST00000471181                   2        41276113               1884
#> 8  ENST00000461221                   1        41276113                 63
#> 9  ENST00000491747                   2        41276113                759
#> 10 ENST00000484087                   1        41256933                498
#> 11 ENST00000478531                   1        41276113                623
#> 12 ENST00000493919                   1        41258543                572
#> 13 ENST00000472490                  NA              NA                 NA
#> 14 ENST00000487825                   1        41256933                266
#> 15 ENST00000461574                   1        41243841                242
#> 16 ENST00000467274                  NA              NA                 NA
#> 17 ENST00000412061                   3        41247883                437
#> 18 ENST00000470026                   1        41276113                649
#> 19 ENST00000477152                   1        41276113                622
#> 20 ENST00000492859                   1        41276113                 59
#> 21 ENST00000497488                   1        41246659                177
#> 22 ENST00000494123                   1        41276113                473
#> 23 ENST00000473961                   1        41256908                319
#> 24 ENST00000476777                   1        41276113                222
#> 25 ENST00000461798                   1        41276113                 63
#> 26 ENST00000489037                   1        41276113                 98
#> 27 ENST00000354071                   6        41276113               1598
#> 28 ENST00000352993                   5        41276113                721
#> 29 ENST00000346315                   4        41276113               1624
#> 30 ENST00000351666                   3        41276113                680
#> 31 ENST00000309486                   4        41246659               1567
#>    Translation.start Translation.species Translation.db_type
#> 1           41197695               human                core
#> 2           41197801               human                core
#> 3           41197695               human                core
#> 4           41197695               human                core
#> 5           41197695               human                core
#> 6           41197695               human                core
#> 7           41197695               human                core
#> 8           41256972               human                core
#> 9           41197695               human                core
#> 10          41215361               human                core
#> 11          41215361               human                core
#> 12          41215377               human                core
#> 13                NA                <NA>                <NA>
#> 14          41228505               human                core
#> 15          41228554               human                core
#> 16                NA                <NA>                <NA>
#> 17          41245587               human                core
#> 18          41245601               human                core
#> 19          41245603               human                core
#> 20          41262552               human                core
#> 21          41246129               human                core
#> 22          41246129               human                core
#> 23          41246187               human                core
#> 24          41247863               human                core
#> 25          41256972               human                core
#> 26          41256206               human                core
#> 27          41197695               human                core
#> 28          41197695               human                core
#> 29          41197695               human                core
#> 30          41197695               human                core
#> 31          41197695               human                core
#>    Translation.object_type  Translation.id Translation.Parent strand
#> 1              Translation ENSP00000350283    ENST00000357654     -1
#> 2              Translation ENSP00000417148    ENST00000468300     -1
#> 3              Translation ENSP00000465818    ENST00000586385     -1
#> 4              Translation ENSP00000467329    ENST00000591534     -1
#> 5              Translation ENSP00000465347    ENST00000591849     -1
#> 6              Translation ENSP00000418775    ENST00000493795     -1
#> 7              Translation ENSP00000418960    ENST00000471181     -1
#> 8              Translation ENSP00000418548    ENST00000461221     -1
#> 9              Translation ENSP00000420705    ENST00000491747     -1
#> 10             Translation ENSP00000419481    ENST00000484087     -1
#> 11             Translation ENSP00000420412    ENST00000478531     -1
#> 12             Translation ENSP00000418819    ENST00000493919     -1
#> 13                    <NA>            <NA>               <NA>     -1
#> 14             Translation ENSP00000418212    ENST00000487825     -1
#> 15             Translation ENSP00000417241    ENST00000461574     -1
#> 16                    <NA>            <NA>               <NA>     -1
#> 17             Translation ENSP00000397145    ENST00000412061     -1
#> 18             Translation ENSP00000419274    ENST00000470026     -1
#> 19             Translation ENSP00000419988    ENST00000477152     -1
#> 20             Translation ENSP00000420253    ENST00000492859     -1
#> 21             Translation ENSP00000418986    ENST00000497488     -1
#> 22             Translation ENSP00000419103    ENST00000494123     -1
#> 23             Translation ENSP00000420201    ENST00000473961     -1
#> 24             Translation ENSP00000417554    ENST00000476777     -1
#> 25             Translation ENSP00000417988    ENST00000461798     -1
#> 26             Translation ENSP00000420781    ENST00000489037     -1
#> 27             Translation ENSP00000326002    ENST00000354071     -1
#> 28             Translation ENSP00000312236    ENST00000352993     -1
#> 29             Translation ENSP00000246907    ENST00000346315     -1
#> 30             Translation ENSP00000338007    ENST00000351666     -1
#> 31             Translation ENSP00000310938    ENST00000309486     -1
#>            Exon is_canonical seq_region_name species         source db_type
#> 1  c("core"....            0              17   human ensembl_havana    core
#> 2  c("17", ....            0              17   human ensembl_havana    core
#> 3  c(1, 1, ....            0              17   human         havana    core
#> 4  c("17", ....            0              17   human         havana    core
#> 5  c("17", ....            0              17   human         havana    core
#> 6  c("ENSE0....            0              17   human ensembl_havana    core
#> 7  c(1, 1, ....            1              17   human ensembl_havana    core
#> 8  c("ENSE0....            0              17   human         havana    core
#> 9  c("core"....            0              17   human         havana    core
#> 10 c("ENSE0....            0              17   human         havana    core
#> 11 c(1, 1, ....            0              17   human         havana    core
#> 12 c("ENSE0....            0              17   human         havana    core
#> 13 c("17", ....            0              17   human         havana    core
#> 14 c("core"....            0              17   human         havana    core
#> 15 c("Exon"....            0              17   human         havana    core
#> 16 c("ENSE0....            0              17   human         havana    core
#> 17 c(-1, -1....            0              17   human         havana    core
#> 18 c(412771....            0              17   human         havana    core
#> 19 c(412773....            0              17   human         havana    core
#> 20 c(1, 1, ....            0              17   human         havana    core
#> 21 c(412773....            0              17   human         havana    core
#> 22 c("ENSE0....            0              17   human         havana    core
#> 23 c("ENSE0....            0              17   human         havana    core
#> 24 c(412772....            0              17   human         havana    core
#> 25 c("ENSE0....            0              17   human         havana    core
#> 26 c("ENSE0....            0              17   human         havana    core
#> 27 c("17", ....            0              17   human        ensembl    core
#> 28 c(412775....            0              17   human        ensembl    core
#> 29 c("core"....            0              17   human        ensembl    core
#> 30 c("ENSE0....            0              17   human        ensembl    core
#> 31 c(1, 1, ....            0              17   human        ensembl    core
#>    object_type          Parent      end length display_name version    start
#> 1   Transcript ENSG00000012048 41277387   7094    BRCA1-001       3 41196312
#> 2   Transcript ENSG00000012048 41277468   3273    BRCA1-007       1 41196822
#> 3   Transcript ENSG00000012048 41277346    781    BRCA1-023       1 41197580
#> 4   Transcript ENSG00000012048 41277346   1282    BRCA1-024       1 41197580
#> 5   Transcript ENSG00000012048 41277346    563    BRCA1-025       1 41197580
#> 6   Transcript ENSG00000012048 41277419   5732    BRCA1-006       1 41197646
#> 7   Transcript ENSG00000012048 41277500   5936    BRCA1-005       2 41197646
#> 8   Transcript ENSG00000012048 41277305   5693    BRCA1-010       1 41197695
#> 9   Transcript ENSG00000012048 41277373   2379    BRCA1-014       2 41197695
#> 10  Transcript ENSG00000012048 41256933   1495    BRCA1-015       1 41215361
#> 11  Transcript ENSG00000012048 41277376   1972    BRCA1-009       1 41215361
#> 12  Transcript ENSG00000012048 41277419   1948    BRCA1-008       1 41215377
#> 13  Transcript ENSG00000012048 41223083    561    BRCA1-021       1 41219291
#> 14  Transcript ENSG00000012048 41256933    800    BRCA1-019       1 41228505
#> 15  Transcript ENSG00000012048 41243841    726    BRCA1-022       1 41228554
#> 16  Transcript ENSG00000012048 41277332   4497    BRCA1-012       1 41243115
#> 17  Transcript ENSG00000012048 41247883   1312    BRCA1-026       3 41245587
#> 18  Transcript ENSG00000012048 41277340   2108    BRCA1-011       1 41245601
#> 19  Transcript ENSG00000012048 41277381   1980    BRCA1-004       1 41245603
#> 20  Transcript ENSG00000012048 41277317   1584    BRCA1-002       1 41246129
#> 21  Transcript ENSG00000012048 41277317    779    BRCA1-003       1 41246129
#> 22  Transcript ENSG00000012048 41277467   1612    BRCA1-013       1 41246129
#> 23  Transcript ENSG00000012048 41256908    958    BRCA1-018       1 41246187
#> 24  Transcript ENSG00000012048 41277370    769    BRCA1-017       1 41247863
#> 25  Transcript ENSG00000012048 41277387    582    BRCA1-020       1 41251848
#> 26  Transcript ENSG00000012048 41277338    455    BRCA1-016       1 41256206
#> 27  Transcript ENSG00000012048 41277500   6411    BRCA1-205       3 41196313
#> 28  Transcript ENSG00000012048 41277500   3780    BRCA1-204       3 41196313
#> 29  Transcript ENSG00000012048 41277468   6451    BRCA1-202       3 41196313
#> 30  Transcript ENSG00000012048 41276132   3444    BRCA1-203       3 41196313
#> 31  Transcript ENSG00000012048 41277468   7114    BRCA1-201       4 41196313
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
#> $seq_region_name
#> [1] "17"
#> 
#> $strand
#> [1] -1
#> 
#> $logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> $biotype
#> [1] "protein_coding"
#> 
#> $id
#> [1] "ENSG00000012048"
#> 
#> $assembly_name
#> [1] "GRCh37"
#> 
#> $start
#> [1] 41196312
#> 
#> $version
#> [1] 15
#> 
#> $end
#> [1] 41277500
#> 
#> $canonical_transcript
#> [1] "ENST00000471181.2"
#> 
#> $display_name
#> [1] "BRCA1"
#> 
#> $description
#> [1] "breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]"
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
```
