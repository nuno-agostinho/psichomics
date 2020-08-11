context("Loading SRA data using recount")

test_that("Load SRA project", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()
    
    dataTypes <- c("Gene expression", "Junction quantification",
                   "Sample metadata")
    
    sraSRP050193 <- loadSRAproject("SRP050193")[[1]]
    expect_is(sraSRP050193, "list")
    expect_equal(names(sraSRP050193), dataTypes)
    expect_is(sraSRP050193$`Gene expression`,         "data.frame")
    expect_is(sraSRP050193$`Junction quantification`, "data.frame")
    expect_is(sraSRP050193$`Sample metadata`,         "data.frame")
    expect_is(sraSRP050193$`Sample metadata`$`cell line`, "character")
    expect_true(all(
        startsWith(sraSRP050193$`Sample metadata`$`cell line`, "SUM")))
    
    # Properly format extra information
    sraSRP042620 <- loadSRAproject("SRP042620")[[1]]
    expect_is(sraSRP042620, "list")
    expect_equal(names(sraSRP042620), dataTypes)
    expect_is(sraSRP042620$`Gene expression`,         "data.frame")
    expect_is(sraSRP042620$`Junction quantification`, "data.frame")
    expect_is(sraSRP042620$`Sample metadata`,         "data.frame")
    expect_is(sraSRP042620$`Sample metadata`$tissue, "character")
    expect_true(all(
        grepl("Breast|cancer", sraSRP042620$`Sample metadata`$tissue)))
})
