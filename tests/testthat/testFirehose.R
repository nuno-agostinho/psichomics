context("Retrieve information from Firehose API and download files")

library(jsonlite)
library(httr)

test_that("isFirehoseUp checks if Firehose API is running", {
    link <- paste0("http://firebrowse.org/api/v1/Metadata/HeartBeat")
    heartbeat <- GET(link, query = list(format = "json"))
    up <- isFirehoseUp()
    
    if (status_code(heartbeat) == 200)
        expect_true(up)
    else
        expect_warning(up)
})

# Only test these functions if the API is running
if (isFirehoseUp()) {
    test_that("queryFirehoseData queries the Firehose API", {
        # Query one cohort
        cohort <- "ACC"
        resp <- queryFirehoseData(cohort = cohort)
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
        resp <- queryFirehoseData(cohort = cohorts)
        expect_is(resp, "response")
        parsed <- fromJSON(content(resp, "text"))$StandardData
        expect_is(parsed, "data.frame")
        expect_equal(unique(parsed$cohort), cohorts)
    })
    
    test_that("parseFirehoseMetadata parses metadata from Firehose", {
        # Parse cohorts metadata
        parsed <- parseFirehoseMetadata("Cohorts")
        expect_is(parsed, "list")
        expect_is(parsed[[1]], "data.frame")
        cohorts <- c("ACC", "BLCA", "COAD", "GBM")
        expect_true(all(cohorts %in% parsed[[1]][, 1]))
        
        # Parse centers metadata
        parsed <- parseFirehoseMetadata("Centers")
        expect_is(parsed, "list")
        expect_is(parsed[[1]], "data.frame")
        centers <- c("genome.wustl.edu", "ucsc.edu", "unc.edu")
        expect_true(all(centers %in% parsed[[1]][, 1]))
 
        # Parse dates metadata
        parsed <- parseFirehoseMetadata("Dates")
        expect_is(parsed, "list")
        expect_is(parsed[[1]], "character")
        dates <- c("2015_11_01", "2015_04_02", "2014_10_17")
        expect_true(all(dates %in% parsed[[1]]))
    })
    
    test_that("getFirehoseDates obtains the datestamps of the data available", {
        resp <- getFirehoseDates()
        expect_is(resp, "character")
        dates <- c("2015_11_01", "2015_04_02", "2014_10_17")
        expect_true(all(dates %in% resp))
    })
    
    test_that("getFirehoseCohorts obtains the cohorts available", {
        # Retrieve all cohorts
        parsed <- getFirehoseCohorts()
        expect_is(parsed, "character")
        expect_true(all(c("Adrenocortical carcinoma" = "ACC",
                          "Brain Lower Grade Glioma" = "LGG") %in% parsed))
        
        # Filter by cohorts of interest
        parsed <- getFirehoseCohorts(c("ACC", "LGG"))
        expect_is(parsed, "character")
        expect_equal(parsed, c("Adrenocortical carcinoma" = "ACC",
                               "Brain Lower Grade Glioma" = "LGG"))
    })
}