# Download and process TCGA data

TCGA data obtained via [FireBrowse](http://firebrowse.org/api-docs/)

## Usage

``` r
loadTCGAdata(
  folder = getDownloadsFolder(),
  data = c("clinical", "junction_quantification", "RSEM_genes"),
  exclude = c(".aux.", ".mage-tab.", "MANIFEST.txt"),
  ...,
  download = TRUE
)
```

## Arguments

- folder:

  Character: directory to store the downloaded archives (by default,
  saves to
  [`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md))

- data:

  Character: data to load (see
  [`getTCGAdataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md))

- exclude:

  Character: files and folders to exclude from downloading and from
  loading into R (by default, exclude files containing `.aux.`,
  `.mage-tab.` and `MANIFEST.TXT`)

- ...:

  Arguments passed on to
  [`queryFirebrowseData`](https://nuno-agostinho.github.io/psichomics/reference/queryFirebrowseData.md)

  `date`

  :   Character: dates of the data retrieval by FireBrowse (by default,
      it uses the most recent data available)

  `cohort`

  :   Character: abbreviation of the cohorts (by default, returns data
      for all cohorts)

  `data_type`

  :   Character: data types (optional)

  `tool`

  :   Character: data produced by the selected FireBrowse tools
      (optional)

  `platform`

  :   Character: data generation platforms (optional)

  `center`

  :   Character: data generation centres (optional)

  `level`

  :   Integer: data levels (optional)

  `protocol`

  :   Character: sample characterization protocols (optional)

  `page`

  :   Integer: page of the results to return (optional)

  `page_size`

  :   Integer: number of records per page of results (optional)

  `sort_by`

  :   String: column used to sort the data (by default, sort by cohort)

- download:

  Boolean: download missing files

## Value

A list with the loaded data, unless required files are unavailable and
`download = FALSE` (if so, it returns the URL of files to download)

## See also

Other functions associated with TCGA data retrieval:
[`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md),
[`getTCGAdataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md),
[`isFirebrowseUp()`](https://nuno-agostinho.github.io/psichomics/reference/isFirebrowseUp.md),
[`parseTCGAsampleTypes()`](https://nuno-agostinho.github.io/psichomics/reference/parseTcgaSampleInfo.md)

Other functions to load data:
[`loadGtexData()`](https://nuno-agostinho.github.io/psichomics/reference/loadGtexData.md),
[`loadLocalFiles()`](https://nuno-agostinho.github.io/psichomics/reference/loadLocalFiles.md),
[`loadSRAproject()`](https://nuno-agostinho.github.io/psichomics/reference/loadSRAproject.md)

## Examples

``` r
getTCGAcohorts()
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
getTCGAdataTypes()
#> $`RNA sequencing`
#>   Junction quantification       Exon quantification           Exon expression 
#> "junction_quantification"     "exon_quantification"         "exon_expression" 
#>       Junction expression                RSEM genes     RSEM genes normalized 
#>     "junction_expression"              "RSEM_genes"   "RSEM_genes_normalized" 
#>             RSEM isoforms                Preprocess 
#>           "RSEM_isoforms"              "Preprocess" 
#> 
if (FALSE) { # \dontrun{
loadTCGAdata(cohort = "ACC", data_type = "Clinical")
} # }
```
