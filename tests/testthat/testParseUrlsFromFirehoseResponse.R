context("Test the parsing of URLs from a Firehose response")

test_that("Parse the URLs from a Firehose response", {
    res <- queryFirehoseData(cohort = "ACC")
    url <- parseUrlsFromFirehoseResponse(res)
    expect_equal(length(unique(names(url))), 1)
    expect_true(grepl("Adrenocortical carcinoma", unique(names(url))))
})