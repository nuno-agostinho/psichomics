context("Filter groups according to the number of non-missing values")

test_that("Groups are ignored with less elements than the threshold", {
    vector <- 1:12
    group <- c("yellow", "yellow", rep(c("red", "blue"), 5))

    # Do not allow yellow
    filtered <- filterGroups(vector, group, threshold=3)
    
    names(vector) <- group
    # Order by "red" than "blue" and ignore "yellow"
    ordered <- vector[3:12][order(names(vector[3:12]), decreasing = TRUE)]
    expect_identical(filtered, ordered)
})

test_that("Groups are maintained with same elements as the threshold", {
    vector <- 1:12
    group <- c("yellow", "yellow", rep(c("red", "blue"), 5))
    
    # Allow all groups
    filtered <- filterGroups(vector, group, threshold=2)
    
    names(vector) <- group
    # Order by "yellow", "red" than "blue"
    ordered <- vector[order(names(vector), decreasing = TRUE)]
    expect_identical(filtered, ordered)
})

test_that("Groups are ignored with less non-missing values than the threshold", {
    vector <- c(NA, 2:12)
    group <- c("yellow", "yellow", rep(c("red", "blue"), 5))
    
    # Do not allow yellow
    filtered <- filterGroups(vector, group, threshold=2)
    
    names(vector) <- group
    # Order by "red" than "blue" and ignore "yellow"
    ordered <- vector[3:12][order(names(vector[3:12]), decreasing = TRUE)]
    expect_identical(filtered, ordered)
})

