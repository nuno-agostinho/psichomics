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
#> [1] "22928"
#> 
#> $search$retmax
#> [1] "3"
#> 
#> $search$retstart
#> [1] "0"
#> 
#> $search$idlist
#> [1] "29687286" "38421676" "29470806"
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
#> [1] "(\"brca1 protein human\"[Supplementary Concept] OR \"brca1 protein human\"[All Fields] OR \"brca1\"[All Fields] OR \"genes brca1\"[MeSH Terms] OR (\"genes\"[All Fields] AND \"brca1\"[All Fields]) OR \"brca1 genes\"[All Fields] OR \"brca1 s\"[All Fields]) AND (\"cancer s\"[All Fields] OR \"cancerated\"[All Fields] OR \"canceration\"[All Fields] OR \"cancerization\"[All Fields] OR \"cancerized\"[All Fields] OR \"cancerous\"[All Fields] OR \"neoplasms\"[MeSH Terms] OR \"neoplasms\"[All Fields] OR \"cancer\"[All Fields] OR \"cancers\"[All Fields] OR (\"adrenocortical carcinoma\"[MeSH Terms] OR (\"adrenocortical\"[All Fields] AND \"carcinoma\"[All Fields]) OR \"adrenocortical carcinoma\"[All Fields]))"
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
#> $`29470806`
#> $`29470806`$uid
#> [1] "29470806"
#> 
#> $`29470806`$pubdate
#> [1] "2018 Jul"
#> 
#> $`29470806`$epubdate
#> [1] "2018 Feb 22"
#> 
#> $`29470806`$source
#> [1] "Breast Cancer Res Treat"
#> 
#> $`29470806`$authors
#>                name authtype clusterid
#> 1           Singh J   Author          
#> 2           Thota N   Author          
#> 3           Singh S   Author          
#> 4           Padhi S   Author          
#> 5           Mohan P   Author          
#> 6         Deshwal S   Author          
#> 7             Sur S   Author          
#> 8           Ghosh M   Author          
#> 9         Agarwal A   Author          
#> 10          Sarin R   Author          
#> 11          Ahmed R   Author          
#> 12          Almel S   Author          
#> 13    Chakraborti B   Author          
#> 14          Raina V   Author          
#> 15     DadiReddy PK   Author          
#> 16        Smruti BK   Author          
#> 17        Rajappa S   Author          
#> 18     Dodagoudar C   Author          
#> 19       Aggarwal S   Author          
#> 20        Singhal M   Author          
#> 21          Joshi A   Author          
#> 22          Kumar R   Author          
#> 23          Kumar A   Author          
#> 24        Mishra DK   Author          
#> 25          Arora N   Author          
#> 26         Karaba A   Author          
#> 27       Sankaran S   Author          
#> 28     Katragadda S   Author          
#> 29          Ghosh A   Author          
#> 30 Veeramachaneni V   Author          
#> 31      Hariharan R   Author          
#> 32        Mannan AU   Author          
#> 
#> $`29470806`$lastauthor
#> [1] "Mannan AU"
#> 
#> $`29470806`$title
#> [1] "Screening of over 1000 Indian patients with breast and/or ovarian cancer with a multi-gene panel: prevalence of BRCA1/2 and non-BRCA mutations."
#> 
#> $`29470806`$sorttitle
#> [1] "screening of over 1000 indian patients with breast and or ovarian cancer with a multi gene panel prevalence of brca1 2 and non brca mutations"
#> 
#> $`29470806`$volume
#> [1] "170"
#> 
#> $`29470806`$issue
#> [1] "1"
#> 
#> $`29470806`$pages
#> [1] "189-196"
#> 
#> $`29470806`$lang
#> [1] "eng"
#> 
#> $`29470806`$nlmuniqueid
#> [1] "8111104"
#> 
#> $`29470806`$issn
#> [1] "0167-6806"
#> 
#> $`29470806`$essn
#> [1] "1573-7217"
#> 
#> $`29470806`$pubtype
#> [1] "Journal Article"
#> 
#> $`29470806`$recordstatus
#> [1] "PubMed - indexed for MEDLINE"
#> 
#> $`29470806`$pubstatus
#> [1] "256"
#> 
#> $`29470806`$articleids
#>   idtype idtypen                     value
#> 1 pubmed       1                  29470806
#> 2    doi       3 10.1007/s10549-018-4726-x
#> 3    pii       4 10.1007/s10549-018-4726-x
#> 
#> $`29470806`$history
#>   pubstatus             date
#> 1  received 2017/11/16 00:00
#> 2  accepted 2018/02/19 00:00
#> 3    pubmed 2018/02/23 06:00
#> 4   medline 2019/03/05 06:00
#> 5    entrez 2018/02/23 06:00
#> 
#> $`29470806`$references
#> list()
#> 
#> $`29470806`$attributes
#> [1] "Has Abstract"
#> 
#> $`29470806`$pmcrefcount
#> [1] ""
#> 
#> $`29470806`$fulljournalname
#> [1] "Breast cancer research and treatment"
#> 
#> $`29470806`$elocationid
#> [1] "doi: 10.1007/s10549-018-4726-x"
#> 
#> $`29470806`$doctype
#> [1] "citation"
#> 
#> $`29470806`$srccontriblist
#> list()
#> 
#> $`29470806`$booktitle
#> [1] ""
#> 
#> $`29470806`$medium
#> [1] ""
#> 
#> $`29470806`$edition
#> [1] ""
#> 
#> $`29470806`$publisherlocation
#> [1] ""
#> 
#> $`29470806`$publishername
#> [1] ""
#> 
#> $`29470806`$srcdate
#> [1] ""
#> 
#> $`29470806`$reportnumber
#> [1] ""
#> 
#> $`29470806`$availablefromurl
#> [1] ""
#> 
#> $`29470806`$locationlabel
#> [1] ""
#> 
#> $`29470806`$doccontriblist
#> list()
#> 
#> $`29470806`$docdate
#> [1] ""
#> 
#> $`29470806`$bookname
#> [1] ""
#> 
#> $`29470806`$chapter
#> [1] ""
#> 
#> $`29470806`$sortpubdate
#> [1] "2018/07/01 00:00"
#> 
#> $`29470806`$sortfirstauthor
#> [1] "Singh J"
#> 
#> $`29470806`$vernaculartitle
#> [1] ""
#> 
#> 
```
