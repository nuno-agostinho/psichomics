context("Retrieve information from Firebrowse web API and download files")

test_that("getGtexDataURL gets URL for GTEx data v8", {
    url <- getGtexDataURL(8)
    expect_length(url, 4)
    expect_equal(names(url), c("sampleInfo", "subjectInfo", "geneExpr", "junctionQuant"))
    expect_equal(file_ext(url), c("txt", "txt", "gz", "gz"))
    expect_match(basename(url), "GTEx_Analysis_")
    expect_match(basename(url), "V8", ignore.case=TRUE)

    expect_is(attr(url, "date"), "Date")
    expect_equal( as.character(attr(url, "date")), "2019-08-26" )
})

test_that("getGtexDataURL gets URL for GTEx data v7", {
    url <- getGtexDataURL(7)
    expect_length(url, 4)
    expect_equal(names(url), c("sampleInfo", "subjectInfo", "geneExpr", "junctionQuant"))
    expect_equal(file_ext(url), c("txt", "txt", "gz", "gz"))
    expect_match(basename(url), "GTEx_")
    expect_match(basename(url), "V7", ignore.case=TRUE)

    expect_is(attr(url, "date"), "Date")
    expect_equal( as.character(attr(url, "date")), "2017-09-05" )
})

test_that("getGtexDataURL gets URL for GTEx data v6", {
    url <- getGtexDataURL(6)
    expect_length(url, 4)
    expect_equal(names(url), c("sampleInfo", "subjectInfo", "geneExpr", "junctionQuant"))
    expect_equal(file_ext(url), c("txt", "txt", "gz", "gz"))
    expect_match(basename(url), "GTEx_")
    expect_match(basename(url), "V6", ignore.case=TRUE)

    expect_is(attr(url, "date"), "Date")
    expect_equal( as.character(attr(url, "date")), "2017-04-25" )
})

test_that("getGtexDataURL gets URL for GTEx data v4", {
    url <- getGtexDataURL(4)
    expect_length(url, 4)
    expect_equal(names(url), c("sampleInfo", "subjectInfo", "geneExpr", "junctionQuant"))
    expect_equal(file_ext(url), c("txt", "txt", "gz", "gz"))
    expect_match(basename(url), "GTEx_")
    expect_match(basename(url), "V4")

    expect_is(attr(url, "date"), "Date")
    expect_equal( as.character(attr(url, "date")), "2017-04-24" )
})
