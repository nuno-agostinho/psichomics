context("Variance by row of matrix")

test_that("Calculate variance for a one-row matrix", {
    mat <- matrix(rnorm(500), nrow = 1)
    res1 <- rowVar(mat)
    res2 <- apply(mat, 1, var) # R way
    expect_identical(signif(res1, 12), signif(res2, 12))
})

test_that("Calculate variance for a one-row matrix ignoring NAs", {
    mat <- matrix(rnorm(500), nrow = 1)
    mat[rnorm(500) > 0] <- NA
    
    # Ignore NAs
    res1 <- rowVar(mat, na.rm = TRUE)
    res2 <- apply(mat, 1, var, na.rm = TRUE) # R way
    expect_identical(signif(res1, 12), signif(res2, 12))
    
    # Don't ignore NAs
    res1 <- rowVar(mat, na.rm = FALSE)
    res2 <- apply(mat, 1, var, na.rm = FALSE) # R way
    expect_identical(signif(res1, 12), signif(res2, 12))
})

test_that("Calculate variance for a multi-row matrix", {
    mat <- replicate(500, rnorm(500))
    res1 <- rowVar(mat)
    res2 <- apply(mat, 1, var) # R way
    expect_identical(signif(res1, 12), signif(res2, 12))
})

test_that("Calculate variance for a multi-row matrix", {
    mat <- replicate(500, rnorm(500))
    mat[replicate(500, rnorm(500)) > 0] <- NA
    
    # Ignore NAs
    res1 <- rowVar(mat, na.rm = TRUE)
    res2 <- apply(mat, 1, var, na.rm = TRUE) # R way
    expect_identical(signif(res1, 12), signif(res2, 12))
    
    # Don't ignore NAs
    res1 <- rowVar(mat, na.rm = FALSE)
    res2 <- apply(mat, 1, var, na.rm = FALSE) # R way
    expect_identical(signif(res1, 12), signif(res2, 12))
})