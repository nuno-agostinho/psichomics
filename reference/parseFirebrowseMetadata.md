# Query the FireBrowse API for metadata

Query the FireBrowse API for metadata

## Usage

``` r
parseFirebrowseMetadata(type, ...)
```

## Arguments

- type:

  Character: metadata to retrieve

- ...:

  Character: parameters to pass to query (optional)

## Value

List with parsed response

## Examples

``` r
psichomics:::parseFirebrowseMetadata("Dates")
#> $Dates
#>  [1] "2016_01_28" "2015_11_01" "2015_08_21" "2015_06_01" "2015_04_02"
#>  [6] "2015_02_04" "2014_12_06" "2014_10_17" "2014_09_02" "2014_07_15"
#> [11] "2014_05_18" "2014_04_16" "2014_03_16"
#> 
psichomics:::parseFirebrowseMetadata("Centers")
#> $Centers
#>                     center center_type code
#> 1           vanderbilt.edu        CGCC   27
#> 2                  jhu.edu        CGCC   28
#> 3                  pnl.gov        CGCC   29
#> 4         genome.wustl.edu        CGCC   30
#> 5                 bcgsc.ca        CGCC   31
#> 6             sanger.ac.uk         GSC   32
#> 7            broad.mit.edu        CGCC   01
#> 8          hms.harvard.edu        CGCC   02
#> 9                  lbl.gov        CGCC   03
#> 10               mskcc.org        CGCC   04
#> 11             jhu-usc.edu        CGCC   05
#> 12         hudsonalpha.org        CGCC   06
#> 13                 unc.edu        CGCC   07
#> 14           broad.mit.edu         GSC   08
#> 15        genome.wustl.edu         GSC   09
#> 16            hgsc.bcm.edu         GSC   10
#> 17     rubicongenomics.com         COM   11
#> 18            hgsc.bcm.edu        CGCC   12
#> 19                bcgsc.ca        CGCC   13
#> 20      broadinstitute.org        GDAC   14
#> 21      systemsbiology.org        GDAC   15
#> 22                 lbl.gov        GDAC   16
#> 23               mskcc.org        GDAC   17
#> 24                ucsc.edu        GDAC   18
#> 25          mdanderson.org        GDAC   19
#> 26          mdanderson.org        CGCC   20
#> 27        genome.wustl.edu        CGCC   21
#> 28              intgen.org        CGCC   22
#> 29 nationwidechildrens.org        CGCC   23
#> 30          mdanderson.org        CGCC   24
#> 31                ucsc.edu         GSC   25
#> 32          mdanderson.org        CGCC   26
#>                                           display_name short_name sort
#> 1                     Vanderbilt University Proteomics       VUMC code
#> 2              The Johns Hopkins University Proteomics        JHU code
#> 3                       Pacific Northwest National Lab       PNNL code
#> 4  Washington University School of Medicine Proteomics       WUSM code
#> 5        Canada's Michael Smith Genome Sciences Centre      BCGSC code
#> 6                      Wellcome Trust Sanger Institute     SANGER code
#> 7                   Broad Institute of MIT and Harvard         BI code
#> 8                               Harvard Medical School        HMS code
#> 9                Lawrence Berkeley National Laboratory        LBL code
#> 10              Memorial Sloan-Kettering Cancer Center      MSKCC code
#> 11   Johns Hopkins / University of Southern California    JHU_USC code
#> 12             HudsonAlpha Institute for Biotechnology       HAIB code
#> 13                        University of North Carolina        UNC code
#> 14                  Broad Institute of MIT and Harvard         BI code
#> 15            Washington University School of Medicine       WUSM code
#> 16                          Baylor College of Medicine        BCM code
#> 17                                    Rubicon Genomics         RG code
#> 18                          Baylor College of Medicine        BCM code
#> 19       Canada's Michael Smith Genome Sciences Centre      BCGSC code
#> 20                  Broad Institute of MIT and Harvard         BI code
#> 21                       Institute for Systems Biology        ISB code
#> 22                Lawrence Berkely National Laboratory        LBL code
#> 23              Memorial Sloan-Kettering Cancer Center      MSKCC code
#> 24                University of California, Santa Cruz       UCSC code
#> 25                                         MD Anderson        MDA code
#> 26       MD Anderson - RPPA Core Facility (Proteomics)        MDA code
#> 27            Washington University School of Medicine       WUSM code
#> 28                                                 IGC        IGC code
#> 29                                             NCH BCR        NCH code
#> 30       MD Anderson - Pathology/Lab Medicine Hamilton        MDA code
#> 31                University of California, Santa Cruz       UCSC code
#> 32  MD Anderson - Institute for Applied Cancer Science        MDA code
#> 
psichomics:::parseFirebrowseMetadata("HeartBeat")
#> $HeartBeat
#> [1] "FireBrowse API at firebrowse.org:8000 is alive"                
#> [2] "Version: 1.1.40 (2019-10-13 13:15:04 c66e6f910b6a89397a4de26c)"
#> [3] "Root Dir: /local/firebrowse/firebrowse_1.1.40"                 
#> [4] "Launched On: 2025_12_06 14:23:53 EST\n"                        
#> 

# Get the abbreviation and description of all cohorts available
psichomics:::parseFirebrowseMetadata("Cohorts")
#> $Cohorts
#>      cohort                                                      description
#> 1       ACC                                         Adrenocortical carcinoma
#> 2      BLCA                                     Bladder Urothelial Carcinoma
#> 3      BRCA                                        Breast invasive carcinoma
#> 4      CESC Cervical squamous cell carcinoma and endocervical adenocarcinoma
#> 5      CHOL                                               Cholangiocarcinoma
#> 6      COAD                                             Colon adenocarcinoma
#> 7  COADREAD                                        Colorectal adenocarcinoma
#> 8      DLBC                  Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
#> 9      ESCA                                            Esophageal carcinoma 
#> 10     FPPP                                              FFPE Pilot Phase II
#> 11      GBM                                          Glioblastoma multiforme
#> 12   GBMLGG                                                           Glioma
#> 13     HNSC                            Head and Neck squamous cell carcinoma
#> 14     KICH                                               Kidney Chromophobe
#> 15    KIPAN                               Pan-kidney cohort (KICH+KIRC+KIRP)
#> 16     KIRC                                Kidney renal clear cell carcinoma
#> 17     KIRP                            Kidney renal papillary cell carcinoma
#> 18     LAML                                           Acute Myeloid Leukemia
#> 19      LGG                                         Brain Lower Grade Glioma
#> 20     LIHC                                   Liver hepatocellular carcinoma
#> 21     LUAD                                              Lung adenocarcinoma
#> 22     LUSC                                     Lung squamous cell carcinoma
#> 23     MESO                                                     Mesothelioma
#> 24       OV                                Ovarian serous cystadenocarcinoma
#> 25     PAAD                                        Pancreatic adenocarcinoma
#> 26     PCPG                               Pheochromocytoma and Paraganglioma
#> 27     PRAD                                          Prostate adenocarcinoma
#> 28     READ                                            Rectum adenocarcinoma
#> 29     SARC                                                          Sarcoma
#> 30     SKCM                                          Skin Cutaneous Melanoma
#> 31     STAD                                           Stomach adenocarcinoma
#> 32     STES                                 Stomach and Esophageal carcinoma
#> 33     TGCT                                      Testicular Germ Cell Tumors
#> 34     THCA                                                Thyroid carcinoma
#> 35     THYM                                                          Thymoma
#> 36     UCEC                             Uterine Corpus Endometrial Carcinoma
#> 37      UCS                                           Uterine Carcinosarcoma
#> 38      UVM                                                   Uveal Melanoma
#> 
# Get the abbreviation and description of the selected cohorts
psichomics:::parseFirebrowseMetadata("Cohorts", cohort = c("ACC", "BRCA"))
#> $Cohorts
#>   cohort               description
#> 1    ACC  Adrenocortical carcinoma
#> 2   BRCA Breast invasive carcinoma
#> 
```
