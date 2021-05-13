context("Filter groups according to the number of non-missing values")

vector <- 1:12
group  <- c("yellow", "yellow", rep(c("red", "blue"), 5))

col2hex <- function(col) rgb(t(col2rgb(col)), maxColorValue=255)
attr(group, "Colour") <- setNames(col2hex(unique(group)), unique(group))

test_that("Groups are ignored with less elements than the threshold", {
    # Discard yellow for having low number of values
    filtered <- filterGroups(vector, group, threshold=3)

    res <- vector[-c(1:2)]

    g <- group[-c(1:2)]
    attr(g, "Colour") <- attr(group, "Colour")
    attr(res, "Groups") <- g
    expect_identical(filtered, res)
})

test_that("Groups are maintained with same elements as the threshold", {
    # Do not discard any of those groups
    filtered <- filterGroups(vector, group, threshold=2)

    res <- vector
    attr(res, "Groups") <- group
    expect_identical(filtered, res)
})

context("Filter groups according to the number of non-missing values")

vector <- c(NA, 2:12)
group  <- c("yellow", "yellow", rep(c("red", "blue"), 5))

test_that("Groups are ignored with less non-missing values than the threshold",{
    # Discard yellow for having low number of non-NA values
    filtered <- filterGroups(vector, group, threshold=2)

    res <- vector[-c(1:2)]
    attr(res, "Groups") <- group[-(1:2)]
    expect_identical(filtered, res)
})
