# Query the PubMed REST API

Query the PubMed REST API

## Usage

``` r
queryPubMed(primary, ..., top = 3, field = "abstract", sort = "relevance")
```

## Arguments

- primary:

  Character: primary search term

- ...:

  Character: other relevant search terms

- top:

  Numeric: number of articles to retrieve

- field:

  Character: field of interest where to look for terms (`abstract` by
  default)

- sort:

  Character: sort by a given parameter (`relevance` by default)

## Value

Parsed response

## Examples

``` r
psichomics:::queryPubMed("BRCA1", "cancer", "adrenocortical carcinoma")
#> $search
#> $search$count
#> [1] "22909"
#> 
#> $search$retmax
#> [1] "3"
#> 
#> $search$retstart
#> [1] "0"
#> 
#> $search$idlist
#> [1] "29687286" "38421676" "33406487"
#> 
#> $search$translationset
#>                       from
#> 1                    BRCA1
#> 2                   cancer
#> 3 adrenocortical carcinoma
#>                                                                                                                                                                                                                                                                             to
#> 1                              "brca1 protein, human"[Supplementary Concept] OR "brca1 protein, human"[All Fields] OR "brca1"[All Fields] OR "genes, brca1"[MeSH Terms] OR ("genes"[All Fields] AND "brca1"[All Fields]) OR "brca1 genes"[All Fields] OR "brca1's"[All Fields]
#> 2 "cancer's"[All Fields] OR "cancerated"[All Fields] OR "canceration"[All Fields] OR "cancerization"[All Fields] OR "cancerized"[All Fields] OR "cancerous"[All Fields] OR "neoplasms"[MeSH Terms] OR "neoplasms"[All Fields] OR "cancer"[All Fields] OR "cancers"[All Fields]
#> 3                                                                                                                               "adrenocortical carcinoma"[MeSH Terms] OR ("adrenocortical"[All Fields] AND "carcinoma"[All Fields]) OR "adrenocortical carcinoma"[All Fields]
#> 
#> $search$querytranslation
#> [1] "(\"brca1 protein human\"[Supplementary Concept] OR \"brca1 protein human\"[All Fields] OR \"brca1\"[All Fields] OR \"genes, brca1\"[MeSH Terms] OR (\"genes\"[All Fields] AND \"brca1\"[All Fields]) OR \"brca1 genes\"[All Fields] OR \"brca1 s\"[All Fields]) AND (\"cancer s\"[All Fields] OR \"cancerated\"[All Fields] OR \"canceration\"[All Fields] OR \"cancerization\"[All Fields] OR \"cancerized\"[All Fields] OR \"cancerous\"[All Fields] OR \"neoplasms\"[MeSH Terms] OR \"neoplasms\"[All Fields] OR \"cancer\"[All Fields] OR \"cancers\"[All Fields] OR (\"adrenocortical carcinoma\"[MeSH Terms] OR (\"adrenocortical\"[All Fields] AND \"carcinoma\"[All Fields]) OR \"adrenocortical carcinoma\"[All Fields]))"
#> 
#> 
#> $`29687286`
#> $`29687286`$uid
#> [1] "29687286"
#> 
#> $`29687286`$pubdate
#> [1] "2020"
#> 
#> $`29687286`$epubdate
#> [1] ""
#> 
#> $`29687286`$source
#> [1] "Adv Exp Med Biol"
#> 
#> $`29687286`$authors
#>         name authtype clusterid
#> 1   Saleem M   Author          
#> 2 Ghazali MB   Author          
#> 3 Wahab MAMA   Author          
#> 4  Yusoff NM   Author          
#> 5   Mahsin H   Author          
#> 6    Seng CE   Author          
#> 7  Khalid IA   Author          
#> 8 Rahman MNG   Author          
#> 9  Yahaya BH   Author          
#> 
#> $`29687286`$lastauthor
#> [1] "Yahaya BH"
#> 
#> $`29687286`$title
#> [1] "The BRCA1 and BRCA2 Genes in Early-Onset Breast Cancer Patients."
#> 
#> $`29687286`$sorttitle
#> [1] "brca1 and brca2 genes in early onset breast cancer patients"
#> 
#> $`29687286`$volume
#> [1] "1292"
#> 
#> $`29687286`$issue
#> [1] ""
#> 
#> $`29687286`$pages
#> [1] "1-12"
#> 
#> $`29687286`$lang
#> [1] "eng"
#> 
#> $`29687286`$nlmuniqueid
#> [1] "0121103"
#> 
#> $`29687286`$issn
#> [1] "0065-2598"
#> 
#> $`29687286`$essn
#> [1] ""
#> 
#> $`29687286`$pubtype
#> [1] "Journal Article" "Review"         
#> 
#> $`29687286`$recordstatus
#> [1] "PubMed - indexed for MEDLINE"
#> 
#> $`29687286`$pubstatus
#> [1] "4"
#> 
#> $`29687286`$articleids
#>   idtype idtypen                 value
#> 1 pubmed       1              29687286
#> 2    doi       3 10.1007/5584_2018_147
#> 
#> $`29687286`$history
#>   pubstatus             date
#> 1    pubmed 2018/04/25 06:00
#> 2   medline 2021/01/26 06:00
#> 3    entrez 2018/04/25 06:00
#> 
#> $`29687286`$references
#> list()
#> 
#> $`29687286`$attributes
#> [1] "Has Abstract"
#> 
#> $`29687286`$pmcrefcount
#> [1] ""
#> 
#> $`29687286`$fulljournalname
#> [1] "Advances in experimental medicine and biology"
#> 
#> $`29687286`$elocationid
#> [1] "doi: 10.1007/5584_2018_147"
#> 
#> $`29687286`$doctype
#> [1] "citation"
#> 
#> $`29687286`$srccontriblist
#> list()
#> 
#> $`29687286`$booktitle
#> [1] ""
#> 
#> $`29687286`$medium
#> [1] ""
#> 
#> $`29687286`$edition
#> [1] ""
#> 
#> $`29687286`$publisherlocation
#> [1] ""
#> 
#> $`29687286`$publishername
#> [1] ""
#> 
#> $`29687286`$srcdate
#> [1] ""
#> 
#> $`29687286`$reportnumber
#> [1] ""
#> 
#> $`29687286`$availablefromurl
#> [1] ""
#> 
#> $`29687286`$locationlabel
#> [1] ""
#> 
#> $`29687286`$doccontriblist
#> list()
#> 
#> $`29687286`$docdate
#> [1] ""
#> 
#> $`29687286`$bookname
#> [1] ""
#> 
#> $`29687286`$chapter
#> [1] ""
#> 
#> $`29687286`$sortpubdate
#> [1] "2020/01/01 00:00"
#> 
#> $`29687286`$sortfirstauthor
#> [1] "Saleem M"
#> 
#> $`29687286`$vernaculartitle
#> [1] ""
#> 
#> 
#> $`38421676`
#> $`38421676`$uid
#> [1] "38421676"
#> 
#> $`38421676`$pubdate
#> [1] "2024 Apr 1"
#> 
#> $`38421676`$epubdate
#> [1] ""
#> 
#> $`38421676`$source
#> [1] "JAMA Oncol"
#> 
#> $`38421676`$authors
#>                                             name       authtype clusterid
#> 1                                     Lubinski J         Author          
#> 2                                  Kotsopoulos J         Author          
#> 3                                       Moller P         Author          
#> 4                                          Pal T         Author          
#> 5                                        Eisen A         Author          
#> 6                                         Peck L         Author          
#> 7                                      Karlan BY         Author          
#> 8                                       Aeilts A         Author          
#> 9                                          Eng C         Author          
#> 10                                   Bordeleau L         Author          
#> 11                                    Foulkes WD         Author          
#> 12                                        Tung N         Author          
#> 13                                      Couch FJ         Author          
#> 14                                     Fruscio R         Author          
#> 15                               Ramon Y Cajal T         Author          
#> 16                                     Singer CF         Author          
#> 17                                  Neuhausen SL         Author          
#> 18                                     Zakalik D         Author          
#> 19                                    Cybulski C         Author          
#> 20                                    Gronwald J         Author          
#> 21                                    Huzarski T         Author          
#> 22                                      Stempa K         Author          
#> 23                                      Dungan J         Author          
#> 24                                   Cullinane C         Author          
#> 25                                    Olopade OI         Author          
#> 26                                    Metcalfe K         Author          
#> 27                                         Sun P         Author          
#> 28                                      Narod SA         Author          
#> 29 Hereditary Breast Cancer Clinical Study Group CollectiveName          
#> 
#> $`38421676`$lastauthor
#> [1] "Narod SA"
#> 
#> $`38421676`$title
#> [1] "MRI Surveillance and Breast Cancer Mortality in Women With BRCA1 and BRCA2 Sequence Variations."
#> 
#> $`38421676`$sorttitle
#> [1] "mri surveillance and breast cancer mortality in women with brca1 and brca2 sequence variations"
#> 
#> $`38421676`$volume
#> [1] "10"
#> 
#> $`38421676`$issue
#> [1] "4"
#> 
#> $`38421676`$pages
#> [1] "493-499"
#> 
#> $`38421676`$lang
#> [1] "eng"
#> 
#> $`38421676`$nlmuniqueid
#> [1] "101652861"
#> 
#> $`38421676`$issn
#> [1] "2374-2437"
#> 
#> $`38421676`$essn
#> [1] "2374-2445"
#> 
#> $`38421676`$pubtype
#> [1] "Journal Article"
#> 
#> $`38421676`$recordstatus
#> [1] "PubMed - indexed for MEDLINE"
#> 
#> $`38421676`$pubstatus
#> [1] "4"
#> 
#> $`38421676`$articleids
#>   idtype idtypen                       value
#> 1 pubmed       1                    38421676
#> 2    pmc       8                 PMC10905376
#> 3  pmcid       5        pmc-id: PMC10905376;
#> 4    doi       3 10.1001/jamaoncol.2023.6944
#> 5    pii       4                     2815702
#> 
#> $`38421676`$history
#>     pubstatus             date
#> 1     medline 2024/04/19 06:43
#> 2      pubmed 2024/02/29 12:42
#> 3      entrez 2024/02/29 11:34
#> 4 pmc-release 2024/02/29 00:00
#> 
#> $`38421676`$references
#>                                                                 refsource
#> 1 JAMA Oncol. 2024 Apr 1;10(4):435-436. doi: 10.1001/jamaoncol.2023.5186.
#>      reftype     pmid note
#> 1 Comment in 38421667     
#> 
#> $`38421676`$attributes
#> [1] "Has Abstract"
#> 
#> $`38421676`$pmcrefcount
#> [1] 13
#> 
#> $`38421676`$fulljournalname
#> [1] "JAMA oncology"
#> 
#> $`38421676`$elocationid
#> [1] "doi: 10.1001/jamaoncol.2023.6944"
#> 
#> $`38421676`$doctype
#> [1] "citation"
#> 
#> $`38421676`$srccontriblist
#> list()
#> 
#> $`38421676`$booktitle
#> [1] ""
#> 
#> $`38421676`$medium
#> [1] ""
#> 
#> $`38421676`$edition
#> [1] ""
#> 
#> $`38421676`$publisherlocation
#> [1] ""
#> 
#> $`38421676`$publishername
#> [1] ""
#> 
#> $`38421676`$srcdate
#> [1] ""
#> 
#> $`38421676`$reportnumber
#> [1] ""
#> 
#> $`38421676`$availablefromurl
#> [1] ""
#> 
#> $`38421676`$locationlabel
#> [1] ""
#> 
#> $`38421676`$doccontriblist
#> list()
#> 
#> $`38421676`$docdate
#> [1] ""
#> 
#> $`38421676`$bookname
#> [1] ""
#> 
#> $`38421676`$chapter
#> [1] ""
#> 
#> $`38421676`$sortpubdate
#> [1] "2024/04/01 00:00"
#> 
#> $`38421676`$sortfirstauthor
#> [1] "Lubinski J"
#> 
#> $`38421676`$vernaculartitle
#> [1] ""
#> 
#> 
#> $`33406487`
#> $`33406487`$uid
#> [1] "33406487"
#> 
#> $`33406487`$pubdate
#> [1] "2021 Jan 6"
#> 
#> $`33406487`$epubdate
#> [1] "2021 Jan 6"
#> 
#> $`33406487`$source
#> [1] "J Natl Compr Canc Netw"
#> 
#> $`33406487`$authors
#>                 name       authtype clusterid
#> 1            Daly MB         Author          
#> 2              Pal T         Author          
#> 3           Berry MP         Author          
#> 4            Buys SS         Author          
#> 5          Dickson P         Author          
#> 6         Domchek SM         Author          
#> 7        Elkhanany A         Author          
#> 8         Friedman S         Author          
#> 9          Goggins M         Author          
#> 10         Hutton ML         Author          
#> 11               CGC CollectiveName          
#> 12         Karlan BY         Author          
#> 13            Khan S         Author          
#> 14           Klein C         Author          
#> 15        Kohlmann W         Author          
#> 16               CGC CollectiveName          
#> 17         Kurian AW         Author          
#> 18         Laronga C         Author          
#> 19         Litton JK         Author          
#> 20            Mak JS         Author          
#> 21              LCGC CollectiveName          
#> 22       Menendez CS         Author          
#> 23       Merajver SD         Author          
#> 24       Norquist BS         Author          
#> 25           Offit K         Author          
#> 26       Pederson HJ         Author          
#> 27          Reiser G         Author          
#> 28               CGC CollectiveName          
#> 29 Senter-Jamieson L         Author          
#> 30               CGC CollectiveName          
#> 31        Shannon KM         Author          
#> 32         Shatsky R         Author          
#> 33     Visvanathan K         Author          
#> 34        Weitzel JN         Author          
#> 35           Wick MJ         Author          
#> 36       Wisinski KB         Author          
#> 37       Yurgelun MB         Author          
#> 38         Darlow SD         Author          
#> 39          Dwyer MA         Author          
#> 
#> $`33406487`$lastauthor
#> [1] "Dwyer MA"
#> 
#> $`33406487`$title
#> [1] "Genetic/Familial High-Risk Assessment: Breast, Ovarian, and Pancreatic, Version 2.2021, NCCN Clinical Practice Guidelines in Oncology."
#> 
#> $`33406487`$sorttitle
#> [1] "genetic familial high risk assessment breast ovarian and pancreatic version 2 2021 nccn clinical practice guidelines in oncology"
#> 
#> $`33406487`$volume
#> [1] "19"
#> 
#> $`33406487`$issue
#> [1] "1"
#> 
#> $`33406487`$pages
#> [1] "77-102"
#> 
#> $`33406487`$lang
#> [1] "eng"
#> 
#> $`33406487`$nlmuniqueid
#> [1] "101162515"
#> 
#> $`33406487`$issn
#> [1] "1540-1405"
#> 
#> $`33406487`$essn
#> [1] "1540-1413"
#> 
#> $`33406487`$pubtype
#> [1] "Journal Article"    "Practice Guideline"
#> 
#> $`33406487`$recordstatus
#> [1] "PubMed - indexed for MEDLINE"
#> 
#> $`33406487`$pubstatus
#> [1] "3"
#> 
#> $`33406487`$articleids
#>   idtype idtypen                   value
#> 1 pubmed       1                33406487
#> 2    doi       3 10.6004/jnccn.2021.0001
#> 3    pii       4            jnccnGLS1901
#> 
#> $`33406487`$history
#>   pubstatus             date
#> 1    entrez 2021/01/06 20:05
#> 2    pubmed 2021/01/07 06:00
#> 3   medline 2021/11/06 06:00
#> 
#> $`33406487`$references
#>                                                                    refsource
#> 1 J Natl Compr Canc Netw. 2022 Feb;20(2):xxvi. doi: 10.6004/jnccn.2022.0014.
#> 2  J Natl Compr Canc Netw. 2022 Feb;20(2):xxv. doi: 10.6004/jnccn.2021.7103.
#>      reftype     pmid note
#> 1 Comment in 35130498     
#> 2 Comment in 35130501     
#> 
#> $`33406487`$attributes
#> [1] "Has Abstract"
#> 
#> $`33406487`$pmcrefcount
#> [1] ""
#> 
#> $`33406487`$fulljournalname
#> [1] "Journal of the National Comprehensive Cancer Network : JNCCN"
#> 
#> $`33406487`$elocationid
#> [1] "doi: 10.6004/jnccn.2021.0001"
#> 
#> $`33406487`$doctype
#> [1] "citation"
#> 
#> $`33406487`$srccontriblist
#> list()
#> 
#> $`33406487`$booktitle
#> [1] ""
#> 
#> $`33406487`$medium
#> [1] ""
#> 
#> $`33406487`$edition
#> [1] ""
#> 
#> $`33406487`$publisherlocation
#> [1] ""
#> 
#> $`33406487`$publishername
#> [1] ""
#> 
#> $`33406487`$srcdate
#> [1] ""
#> 
#> $`33406487`$reportnumber
#> [1] ""
#> 
#> $`33406487`$availablefromurl
#> [1] ""
#> 
#> $`33406487`$locationlabel
#> [1] ""
#> 
#> $`33406487`$doccontriblist
#> list()
#> 
#> $`33406487`$docdate
#> [1] ""
#> 
#> $`33406487`$bookname
#> [1] ""
#> 
#> $`33406487`$chapter
#> [1] ""
#> 
#> $`33406487`$sortpubdate
#> [1] "2021/01/06 00:00"
#> 
#> $`33406487`$sortfirstauthor
#> [1] "Daly MB"
#> 
#> $`33406487`$vernaculartitle
#> [1] ""
#> 
#> 
```
