context("Variance by row of matrix")

test_that("customRowVars calculates the variance per row of a matrix", {
    # Matrix input
    mat  <- replicate(10, rnorm(20))
    res  <- apply(mat, 1, var) # R way
    
    res1 <- customRowVars(mat)
    expect_true(all(res1 - res < 10e-16))
    
    res2 <- customRowVars(mat, fast = TRUE)
    expect_true(all(res2 - res < 10e-16))
    
    # Single vector input
    mat  <- mat[1, ]
    res  <- var(mat)
    res1 <- customRowVars(mat)
    expect_equal(res, res1)
    
    res  <- var(mat)
    res2 <- customRowVars(mat, fast = TRUE)
    expect_equal(res, res2)
})

test_that("Calculate variance for a one-row matrix", {
    mat  <- matrix(rnorm(500), nrow = 1)
    res  <- apply(mat, 1, var) # R way
    
    res1 <- customRowVars(mat)
    expect_equal(signif(res, 16), signif(res1, 16))
    
    res2 <- customRowVars(mat, fast = TRUE)
    expect_equal(signif(res, 16), signif(res2, 16))
})

test_that("Calculate variance for a one-row matrix ignoring NAs", {
    mat <- matrix(rnorm(500), nrow = 1)
    mat[rnorm(500) > 0] <- NA
    
    # Ignore NAs
    res  <- apply(mat, 1, var, na.rm = TRUE) # R way
    res1 <- customRowVars(mat, na.rm = TRUE)
    expect_equal(signif(res, 16), signif(res1, 16))
    
    res2 <- customRowVars(mat, na.rm = TRUE, fast = TRUE)
    expect_equal(signif(res, 16), signif(res2, 16))
    
    # Don't ignore NAs
    res  <- apply(mat, 1, var, na.rm = FALSE) # R way
    res1 <- customRowVars(mat, na.rm = FALSE)
    expect_equal(signif(res, 16), signif(res1, 16))
    
    res2 <- customRowVars(mat, na.rm = FALSE, fast = TRUE)
    expect_equal(signif(res, 16), signif(res2, 16))
})

test_that("Calculate variance for a multi-row matrix", {
    mat  <- replicate(500, rnorm(500))
    res  <- apply(mat, 1, var) # R way
    res1 <- customRowVars(mat)
    expect_equal(signif(res, 16), signif(res1, 16))
    
    res2 <- customRowVars(mat, fast = TRUE)
    expect_equal(signif(res, 16), signif(res2, 16))
})

test_that("Calculate variance for a multi-row matrix", {
    mat <- replicate(500, rnorm(500))
    mat[replicate(500, rnorm(500)) > 0] <- NA
    
    # Ignore NAs
    res  <- apply(mat, 1, var, na.rm = TRUE) # R way
    res1 <- customRowVars(mat, na.rm = TRUE)
    expect_equal(signif(res, 16), signif(res1, 16))
    
    res2 <- customRowVars(mat, na.rm = TRUE, fast = TRUE)
    expect_equal(signif(res, 16), signif(res2, 16))
    
    # Don't ignore NAs
    res  <- apply(mat, 1, var, na.rm = FALSE) # R way
    res1 <- customRowVars(mat, na.rm = FALSE)
    expect_equal(signif(res, 16), signif(res1, 16))
    
    res2 <- customRowVars(mat, na.rm = FALSE, fast = TRUE)
    expect_equal(signif(res, 16), signif(res2, 16))
})
