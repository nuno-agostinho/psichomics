context("Rename duplicated values")

test_that("Return unique values in case there are duplicates", {
    check <- c("cat", "mouse", "snake")
    comp <- c("snake", "spider", "dog")
    rep <- renameDuplicated(check, comp)
    expect_equal(rep, c(comp, check[1:2], "snake (1)"))
})

test_that("Return unique values and count + 1 if needed", {
    check <- c("cat", "mouse", "snake")
    comp <- c("snake", "spider", "dog", "snake (2)", "snake (3)")
    rep <- renameDuplicated(check, comp)
    expect_equal(rep, c(comp, check[1:2], "snake (4)"))
})

test_that("Return same values in case no values are duplicated", {
    check <- c("cat", "mouse")
    comp <- c("snake", "spider", "dog")
    rep <- renameDuplicated(check, comp)
    expect_equal(rep, c(comp, check))
})

test_that("Return comp values only if argument check is of length 0", {
    check <- character(0)
    comp <- c("cat", "mouse")
    rep <- renameDuplicated(check, comp)
    expect_equal(rep, comp)
})

test_that("Return check values only if argument comp is of length 0", {
    check <- c("cat", "mouse")
    comp <- character(0)
    expect_error(renameDuplicated(check, comp))
})