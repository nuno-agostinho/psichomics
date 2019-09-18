context("Variance by row of matrix")

test_that("rowVars calculates the variance per row of a matrix", {
    # Passing a matrix
    mat <- replicate(10, rnorm(20))
    res1 <- rowVars(mat)
    res2 <- apply(mat, 1, var) # R way
    expect_true(all(res1 - res2 < 10e-16))
    
    # Passing a single vector
    mat <- mat[1, ]
    expect_equal(rowVars(mat), var(mat))
})

test_that("Calculate variance for a one-row matrix", {
    mat <- matrix(rnorm(500), nrow = 1)
    res1 <- rowVars(mat)
    res2 <- apply(mat, 1, var) # R way
    expect_equal(signif(res1, 16), signif(res2, 16))
})

test_that("Calculate variance for a one-row matrix ignoring NAs", {
    mat <- matrix(rnorm(500), nrow = 1)
    mat[rnorm(500) > 0] <- NA
    
    # Ignore NAs
    res1 <- rowVars(mat, na.rm = TRUE)
    res2 <- apply(mat, 1, var, na.rm = TRUE) # R way
    expect_equal(signif(res1, 16), signif(res2, 16))
    
    # Don't ignore NAs
    res1 <- rowVars(mat, na.rm = FALSE)
    res2 <- apply(mat, 1, var, na.rm = FALSE) # R way
    expect_equal(signif(res1, 16), signif(res2, 16))
})

test_that("Calculate variance for a multi-row matrix", {
    mat <- replicate(500, rnorm(500))
    res1 <- rowVars(mat)
    res2 <- apply(mat, 1, var) # R way
    expect_equal(signif(res1, 16), signif(res2, 16))
})

test_that("Calculate variance for a multi-row matrix", {
    mat <- replicate(500, rnorm(500))
    mat[replicate(500, rnorm(500)) > 0] <- NA
    
    # Ignore NAs
    res1 <- rowVars(mat, na.rm = TRUE)
    res2 <- apply(mat, 1, var, na.rm = TRUE) # R way
    expect_equal(signif(res1, 16), signif(res2, 16))
    
    # Don't ignore NAs
    res1 <- rowVars(mat, na.rm = FALSE)
    res2 <- apply(mat, 1, var, na.rm = FALSE) # R way
    expect_equal(signif(res1, 16), signif(res2, 16))
})
