context("Retrieve information from Firebrowse web API and download files")

library(jsonlite)
library(httr)

skipIfFirebrowseIsUnavailable <- function() {
    if (!isFirebrowseUp())
        skip("API not available")
}

test_that("isFirebrowseUp checks if Firebrowse web API is running", {
    link <- paste0("http://firebrowse.org/api/v1/Metadata/HeartBeat")
    heartbeat <- tryCatch(GET(link, query = list(format="json")), 
                          error=return)
    if ("error" %in% class(heartbeat))
        skip("Couldn't resolve host name")
    
    up <- isFirebrowseUp()
    
    if (status_code(heartbeat) == 200)
        expect_true(up)
    else
        expect_warning(up)
})

# Only test these functions if the API is running
test_that("queryFirebrowseData queries the Firebrowse web API", {
    skipIfFirebrowseIsUnavailable()
    # Query one cohort
    cohort <- "ACC"
    resp <- queryFirebrowseData(cohort = cohort)
    parsed <- fromJSON(content(resp, "text"))
    
    expect_is(resp, "response")
    expect_is(parsed, "list")
    expect_true("StandardData" %in% names(parsed))
    
    parsed <- parsed$StandardData
    expect_is(parsed, "data.frame")
    expect_true("cohort" %in% names(parsed))
    expect_equal(unique(parsed$cohort), cohort)
    
    # Query multiple cohorts
    cohorts <- c("ACC", "BRCA")
    resp <- queryFirebrowseData(cohort = cohorts)
    expect_is(resp, "response")
    parsed <- fromJSON(content(resp, "text"))$StandardData
    expect_is(parsed, "data.frame")
    expect_equal(unique(parsed$cohort), cohorts)
})

test_that("parseFirebrowseMetadata parses metadata from Firebrowse", {
    skipIfFirebrowseIsUnavailable()
    # Parse cohorts metadata
    parsed <- parseFirebrowseMetadata("Cohorts")
    expect_is(parsed, "list")
    expect_is(parsed[[1]], "data.frame")
    cohorts <- c("ACC", "BLCA", "COAD", "GBM")
    expect_true(all(cohorts %in% parsed[[1]][, 1]))
    
    # Parse centers metadata
    parsed <- parseFirebrowseMetadata("Centers")
    expect_is(parsed, "list")
    expect_is(parsed[[1]], "data.frame")
    centers <- c("genome.wustl.edu", "ucsc.edu", "unc.edu")
    expect_true(all(centers %in% parsed[[1]][, 1]))
    
    # Parse dates metadata
    parsed <- parseFirebrowseMetadata("Dates")
    expect_is(parsed, "list")
    expect_is(parsed[[1]], "character")
    dates <- c("2015_11_01", "2015_04_02", "2014_10_17")
    expect_true(all(dates %in% parsed[[1]]))
})

test_that("getFirebrowseDates obtains the datestamps of the data available", {
    skipIfFirebrowseIsUnavailable()
    resp <- getFirebrowseDates()
    expect_is(resp, "Date")
    dates <- c("2015-11-01", "2015-04-02", "2014-10-17")
    expect_true(all(dates %in% as.character(resp)))
})

test_that("getFirebrowseCohorts obtains the cohorts available", {
    skipIfFirebrowseIsUnavailable()
    # Retrieve all cohorts
    parsed <- getFirebrowseCohorts()
    expect_is(parsed, "character")
    expect_true(all(c("ACC" = "Adrenocortical carcinoma",
                      "LGG" = "Brain Lower Grade Glioma") %in% parsed))
    
    # Filter by cohorts of interest
    parsed <- getFirebrowseCohorts(c("ACC", "LGG"))
    expect_is(parsed, "character")
    expect_equal(parsed, c("ACC" = "Adrenocortical carcinoma",
                           "LGG" = "Brain Lower Grade Glioma"))
})

test_that("getFirebrowseDataTypes obtains data types available in Firebrowse", {
    skipIfFirebrowseIsUnavailable()
    types <- getFirebrowseDataTypes()
    expect_true(all(c("exon_expression",
                      "RSEM_genes_normalized",
                      "junction_quantification") %in% types[[1]]))
})

test_that("Parse the URLs from a Firebrowse response", {
    res <- tryCatch(queryFirebrowseData(cohort = "ACC"), error=return)
    if (is(res, "error")) {
        skip("Could not resolve host name")
    } else {
        url <- parseUrlsFromFirebrowseResponse(res)
        expect_equal(length(unique(names(url))), 1)
        expect_true(grepl("Adrenocortical carcinoma", unique(names(url))))
    }
})

# test_that("prepareFirebrowseArchives prepares archives to be loaded", {
#     skipIfFirebrowseIsUnavailable()
#     folder <- "~/Downloads"
#     f <- "gdac.broadinstitute.org_ACC.Merge_Clinical.Level_1.2015110100.0.0"
#     url <- paste0("http://gdac.broadinstitute.org/runs/stddata__2015_11_01",
#                   "/data/ACC/20151101/", f, ".tar.gz", c("", ".md5"))
#     res <- prepareFirebrowseArchives(url, folder, quiet = TRUE)
#     expect_true(res)
#     # Check if prepared folder exists at indicated destination
#     file <- file.path(folder, f)
#     expect_true(file.exists(file))
#     # Remove folder after testing
#     unlink(file, recursive = TRUE)
# })

