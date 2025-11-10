# Get available parameters for TCGA data

Parameters obtained via [FireBrowse](http://firebrowse.org/api-docs/)

## Usage

``` r
getTCGAdataTypes()

getTCGAdates()

getTCGAcohorts(cohort = NULL)
```

## Arguments

- cohort:

  Character: filter results by cohorts (optional)

## Value

Parsed response

## See also

Other functions associated with TCGA data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`isFirebrowseUp()`](https://nuno-agostinho.github.io/psichomics/reference/isFirebrowseUp.md),
[`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md),
[`parseTCGAsampleTypes()`](https://nuno-agostinho.github.io/psichomics/reference/parseTcgaSampleInfo.md)

## Examples

``` r
getTCGAdataTypes()
#> $`RNA sequencing`
#>   Junction quantification       Exon quantification           Exon expression 
#> "junction_quantification"     "exon_quantification"         "exon_expression" 
#>       Junction expression                RSEM genes     RSEM genes normalized 
#>     "junction_expression"              "RSEM_genes"   "RSEM_genes_normalized" 
#>             RSEM isoforms                Preprocess 
#>           "RSEM_isoforms"              "Preprocess" 
#> 
if (isFirebrowseUp()) getTCGAdates()
#>  [1] "2016-01-28" "2015-11-01" "2015-08-21" "2015-06-01" "2015-04-02"
#>  [6] "2015-02-04" "2014-12-06" "2014-10-17" "2014-09-02" "2014-07-15"
#> [11] "2014-05-18" "2014-04-16" "2014-03-16"
if (isFirebrowseUp()) getTCGAcohorts()
#>                                                                ACC 
#>                                         "Adrenocortical carcinoma" 
#>                                                               BLCA 
#>                                     "Bladder Urothelial Carcinoma" 
#>                                                               BRCA 
#>                                        "Breast invasive carcinoma" 
#>                                                               CESC 
#> "Cervical squamous cell carcinoma and endocervical adenocarcinoma" 
#>                                                               CHOL 
#>                                               "Cholangiocarcinoma" 
#>                                                               COAD 
#>                                             "Colon adenocarcinoma" 
#>                                                           COADREAD 
#>                                        "Colorectal adenocarcinoma" 
#>                                                               DLBC 
#>                  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" 
#>                                                               ESCA 
#>                                            "Esophageal carcinoma " 
#>                                                               FPPP 
#>                                              "FFPE Pilot Phase II" 
#>                                                                GBM 
#>                                          "Glioblastoma multiforme" 
#>                                                             GBMLGG 
#>                                                           "Glioma" 
#>                                                               HNSC 
#>                            "Head and Neck squamous cell carcinoma" 
#>                                                               KICH 
#>                                               "Kidney Chromophobe" 
#>                                                              KIPAN 
#>                               "Pan-kidney cohort (KICH+KIRC+KIRP)" 
#>                                                               KIRC 
#>                                "Kidney renal clear cell carcinoma" 
#>                                                               KIRP 
#>                            "Kidney renal papillary cell carcinoma" 
#>                                                               LAML 
#>                                           "Acute Myeloid Leukemia" 
#>                                                                LGG 
#>                                         "Brain Lower Grade Glioma" 
#>                                                               LIHC 
#>                                   "Liver hepatocellular carcinoma" 
#>                                                               LUAD 
#>                                              "Lung adenocarcinoma" 
#>                                                               LUSC 
#>                                     "Lung squamous cell carcinoma" 
#>                                                               MESO 
#>                                                     "Mesothelioma" 
#>                                                                 OV 
#>                                "Ovarian serous cystadenocarcinoma" 
#>                                                               PAAD 
#>                                        "Pancreatic adenocarcinoma" 
#>                                                               PCPG 
#>                               "Pheochromocytoma and Paraganglioma" 
#>                                                               PRAD 
#>                                          "Prostate adenocarcinoma" 
#>                                                               READ 
#>                                            "Rectum adenocarcinoma" 
#>                                                               SARC 
#>                                                          "Sarcoma" 
#>                                                               SKCM 
#>                                          "Skin Cutaneous Melanoma" 
#>                                                               STAD 
#>                                           "Stomach adenocarcinoma" 
#>                                                               STES 
#>                                 "Stomach and Esophageal carcinoma" 
#>                                                               TGCT 
#>                                      "Testicular Germ Cell Tumors" 
#>                                                               THCA 
#>                                                "Thyroid carcinoma" 
#>                                                               THYM 
#>                                                          "Thymoma" 
#>                                                               UCEC 
#>                             "Uterine Corpus Endometrial Carcinoma" 
#>                                                                UCS 
#>                                           "Uterine Carcinosarcoma" 
#>                                                                UVM 
#>                                                   "Uveal Melanoma" 
```
