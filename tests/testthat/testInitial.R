context("Test general functions")

test_that("trimWhitespace does nothing when there's no need to trim", {
    word <- "this is a test"
    expect_equal(trimWhitespace(word), word)
})

test_that("trimWhitespace trims whitespace from a character vector", {
    word <- c("     this         is    a     test         ",
              "another     simple test     right here  ",
              "one               final                 test,          yay")
    res <- c("this is a test", 
             "another simple test right here",
             "one final test, yay")
    expect_equal(trimWhitespace(word), res)
})

test_that("rm.null removes NULL elements from a vector or a list", {
    v1 <- c(1:6, NULL, 2, NULL, 6, 9)
    v2 <- rm.null(v1)
    expect_equal(v2, c(1:6, 2, 6, 9))
    
    l1 <- list(1:3, 6, NULL, 1, NULL, 4:8, NULL)
    l2 <- rm.null(l1)
    expect_equal(l2, list(1:3, 6, 1, 4:8))
})

test_that("rm.null returns the input if there are no NULL elements", {
    v1 <- c(1:6, 2, 6, 9)
    v2 <- rm.null(v1)
    expect_equal(v2, v1)
    
    l1 <- list(1:3, 6, 1, 4:8)
    l2 <- rm.null(l1)
    expect_equal(l2, l1)
})

test_that("rm.null returns NULL for a vector with only NULL elements", {
    v1 <- c(NULL, NULL, NULL)
    v2 <- rm.null(v1)
    expect_equal(v2, NULL)
})

test_that("rm.null returns an empty list for a list with only NULL elements", {
    l1 <- list(NULL, NULL, NULL)
    l2 <- rm.null(l1)
    expect_equal(l2, list())
})