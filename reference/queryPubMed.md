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
#> [1] "22874"
#> 
#> $search$retmax
#> [1] "3"
#> 
#> $search$retstart
#> [1] "0"
#> 
#> $search$idlist
#> [1] "25329591" "29687286" "15546503"
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
#> $`25329591`
#> $`25329591`$uid
#> [1] "25329591"
#> 
#> $`25329591`$pubdate
#> [1] "2015"
#> 
#> $`25329591`$epubdate
#> [1] ""
#> 
#> $`25329591`$source
#> [1] "Anticancer Agents Med Chem"
#> 
#> $`25329591`$authors
#>           name authtype clusterid
#> 1 Romagnolo AP   Author          
#> 2 Romagnolo DF   Author          
#> 3    Selmin OI   Author          
#> 
#> $`25329591`$lastauthor
#> [1] "Selmin OI"
#> 
#> $`25329591`$title
#> [1] "BRCA1 as target for breast cancer prevention and therapy."
#> 
#> $`25329591`$sorttitle
#> [1] "brca1 as target for breast cancer prevention and therapy"
#> 
#> $`25329591`$volume
#> [1] "15"
#> 
#> $`25329591`$issue
#> [1] "1"
#> 
#> $`25329591`$pages
#> [1] "4-14"
#> 
#> $`25329591`$lang
#> [1] "eng"
#> 
#> $`25329591`$nlmuniqueid
#> [1] "101265649"
#> 
#> $`25329591`$issn
#> [1] "1871-5206"
#> 
#> $`25329591`$essn
#> [1] "1875-5992"
#> 
#> $`25329591`$pubtype
#> [1] "Journal Article" "Review"         
#> 
#> $`25329591`$recordstatus
#> [1] "PubMed - indexed for MEDLINE"
#> 
#> $`25329591`$pubstatus
#> [1] "4"
#> 
#> $`25329591`$articleids
#>   idtype idtypen                             value
#> 1 pubmed       1                          25329591
#> 2    doi       3 10.2174/1871520614666141020153543
#> 3    pii       4                  ACAMC-EPUB-62932
#> 
#> $`25329591`$history
#>   pubstatus             date
#> 1  received 2014/09/05 00:00
#> 2   revised 2014/10/14 00:00
#> 3  accepted 2014/10/15 00:00
#> 4    entrez 2014/10/21 06:00
#> 5    pubmed 2014/10/21 06:00
#> 6   medline 2015/08/13 06:00
#> 
#> $`25329591`$references
#> list()
#> 
#> $`25329591`$attributes
#> [1] "Has Abstract"
#> 
#> $`25329591`$pmcrefcount
#> [1] ""
#> 
#> $`25329591`$fulljournalname
#> [1] "Anti-cancer agents in medicinal chemistry"
#> 
#> $`25329591`$elocationid
#> [1] ""
#> 
#> $`25329591`$doctype
#> [1] "citation"
#> 
#> $`25329591`$srccontriblist
#> list()
#> 
#> $`25329591`$booktitle
#> [1] ""
#> 
#> $`25329591`$medium
#> [1] ""
#> 
#> $`25329591`$edition
#> [1] ""
#> 
#> $`25329591`$publisherlocation
#> [1] ""
#> 
#> $`25329591`$publishername
#> [1] ""
#> 
#> $`25329591`$srcdate
#> [1] ""
#> 
#> $`25329591`$reportnumber
#> [1] ""
#> 
#> $`25329591`$availablefromurl
#> [1] ""
#> 
#> $`25329591`$locationlabel
#> [1] ""
#> 
#> $`25329591`$doccontriblist
#> list()
#> 
#> $`25329591`$docdate
#> [1] ""
#> 
#> $`25329591`$bookname
#> [1] ""
#> 
#> $`25329591`$chapter
#> [1] ""
#> 
#> $`25329591`$sortpubdate
#> [1] "2015/01/01 00:00"
#> 
#> $`25329591`$sortfirstauthor
#> [1] "Romagnolo AP"
#> 
#> $`25329591`$vernaculartitle
#> [1] ""
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
#> $`15546503`
#> $`15546503`$uid
#> [1] "15546503"
#> 
#> $`15546503`$pubdate
#> [1] "2004 Nov"
#> 
#> $`15546503`$epubdate
#> [1] ""
#> 
#> $`15546503`$source
#> [1] "Cancer Sci"
#> 
#> $`15546503`$authors
#>        name authtype clusterid
#> 1 Yoshida K   Author          
#> 2    Miki Y   Author          
#> 
#> $`15546503`$lastauthor
#> [1] "Miki Y"
#> 
#> $`15546503`$title
#> [1] "Role of BRCA1 and BRCA2 as regulators of DNA repair, transcription, and cell cycle in response to DNA damage."
#> 
#> $`15546503`$sorttitle
#> [1] "role of brca1 and brca2 as regulators of dna repair transcription and cell cycle in response to dna damage"
#> 
#> $`15546503`$volume
#> [1] "95"
#> 
#> $`15546503`$issue
#> [1] "11"
#> 
#> $`15546503`$pages
#> [1] "866-71"
#> 
#> $`15546503`$lang
#> [1] "eng"
#> 
#> $`15546503`$nlmuniqueid
#> [1] "101168776"
#> 
#> $`15546503`$issn
#> [1] "1347-9032"
#> 
#> $`15546503`$essn
#> [1] "1349-7006"
#> 
#> $`15546503`$pubtype
#> [1] "Journal Article" "Review"         
#> 
#> $`15546503`$recordstatus
#> [1] "PubMed - indexed for MEDLINE"
#> 
#> $`15546503`$pubstatus
#> [1] "4"
#> 
#> $`15546503`$articleids
#>   idtype idtypen                              value
#> 1 pubmed       1                           15546503
#> 2    pmc       8                        PMC11159131
#> 3  pmcid       5               pmc-id: PMC11159131;
#> 4    doi       3 10.1111/j.1349-7006.2004.tb02195.x
#> 
#> $`15546503`$history
#>     pubstatus             date
#> 1      pubmed 2004/11/18 09:00
#> 2     medline 2005/01/15 09:00
#> 3      entrez 2004/11/18 09:00
#> 4 pmc-release 2005/08/19 00:00
#> 
#> $`15546503`$references
#> list()
#> 
#> $`15546503`$attributes
#> [1] "Has Abstract"
#> 
#> $`15546503`$pmcrefcount
#> [1] 62
#> 
#> $`15546503`$fulljournalname
#> [1] "Cancer science"
#> 
#> $`15546503`$elocationid
#> [1] ""
#> 
#> $`15546503`$doctype
#> [1] "citation"
#> 
#> $`15546503`$srccontriblist
#> list()
#> 
#> $`15546503`$booktitle
#> [1] ""
#> 
#> $`15546503`$medium
#> [1] ""
#> 
#> $`15546503`$edition
#> [1] ""
#> 
#> $`15546503`$publisherlocation
#> [1] ""
#> 
#> $`15546503`$publishername
#> [1] ""
#> 
#> $`15546503`$srcdate
#> [1] ""
#> 
#> $`15546503`$reportnumber
#> [1] ""
#> 
#> $`15546503`$availablefromurl
#> [1] ""
#> 
#> $`15546503`$locationlabel
#> [1] ""
#> 
#> $`15546503`$doccontriblist
#> list()
#> 
#> $`15546503`$docdate
#> [1] ""
#> 
#> $`15546503`$bookname
#> [1] ""
#> 
#> $`15546503`$chapter
#> [1] ""
#> 
#> $`15546503`$sortpubdate
#> [1] "2004/11/01 00:00"
#> 
#> $`15546503`$sortfirstauthor
#> [1] "Yoshida K"
#> 
#> $`15546503`$vernaculartitle
#> [1] ""
#> 
#> 
```
