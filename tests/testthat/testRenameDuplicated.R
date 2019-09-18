context("Rename duplicated values")

test_that("Return unique values in case there are duplicates", {
    check <- c("cat", "mouse", "snake")
    comp <- c("snake", "spider", "dog")
    res <- renameDuplicated(check, comp)
    expect_equal(res, c(check[1:2], "snake (1)"))
})

test_that("Return unique values and count + 1 if needed", {
    check <- c("cat", "mouse", "snake")
    comp <- c("snake", "spider", "dog", "snake (2)", "snake (3)")
    res <- renameDuplicated(check, comp)
    expect_equal(res, c(check[1:2], "snake (4)"))
})

test_that("Return same values in case no values are duplicated", {
    check <- c("cat", "mouse")
    comp <- c("snake", "spider", "dog")
    res <- renameDuplicated(check, comp)
    expect_equal(res, c(check))
})

test_that("Return no check values if argument check has length 0", {
    check <- character(0)
    comp <- c("cat", "mouse")
    res <- renameDuplicated(check, comp)
    expect_equal(res, check)
})

test_that("Return the same check values if argument comp is of length 0", {
    check <- c("cat", "mouse")
    comp <- character(0)
    res <- renameDuplicated(check, comp)
    expect_equal(res, check)
    
    check <- c("cat", "mouse", "cat")
    comp <- character(0)
    res <- renameDuplicated(check, comp)
    expect_equal(res, c("cat", "mouse", "cat (1)"))
})

test_that("Return renamed vector in the original order", {
    check <- c("cat", "mouse")
    comp  <- c("cat")
    res <- renameDuplicated(check, comp)
    expect_equal(res, c("cat (1)", "mouse"))
})

test_that("Return renamed vector with no duplicates", {
    check <- c("cat", "mouse", "cat")
    comp  <- c("cat")
    res <- renameDuplicated(check, comp)
    expect_equal(res, c("cat (1)", "mouse", "cat (2)"))
})
