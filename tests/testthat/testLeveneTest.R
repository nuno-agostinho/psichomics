context("Perform Levene's test")

data <- mtcars
values <- data$mpg
groups <- data$carb

test_that("Calculate spread using the median values per group", {
    lev <- leveneTest(values, groups)
    car <- car::leveneTest(values, groups)
    
    expect_is(lev, "htest")
    expect_equal(lev$statistic[[1]], car$`F value`[[1]])
    expect_equal(lev$p.value, car$`Pr(>F)`[[1]])
    expect_true(grepl("median", lev$method))
})

test_that("Calculate spread using the mean values per group", {
    lev <- leveneTest(values, groups, "mean")
    car <- car::leveneTest(values, groups, "mean")
    
    expect_is(lev, "htest")
    expect_equal(lev$statistic[[1]], car$`F value`[[1]])
    expect_equal(lev$p.value, car$`Pr(>F)`[[1]])
    expect_true(grepl("mean", lev$method))
})


test_that("Remove missing values", {
    random <- round(runif(5, 1, length(values)))
    values[random] <- NA
    
    lev <- leveneTest(values, groups)
    car <- car::leveneTest(values, groups)
    
    expect_is(lev, "htest")
    expect_equal(lev$statistic[[1]], car$`F value`[[1]])
    expect_equal(lev$p.value, car$`Pr(>F)`[[1]])
    expect_true(grepl("median", lev$method))
})

test_that("Factorise groups", {
    lev_factor <- leveneTest(values, groups)
    expect_is(lev_factor, "htest")
    
    groups_char <- as.character(groups)
    
    lev_char <- leveneTest(values, groups_char)
    expect_is(lev_char, "htest")
    expect_equal(lev_factor$statistic, lev_char$statistic)
    expect_equal(lev_factor$p.value, lev_char$p.value)
    
    car <- suppressWarnings( car::leveneTest(values, groups_char) )
    expect_equal(lev_char$statistic[[1]], car$`F value`[[1]])
    expect_equal(lev_char$p.value, car$`Pr(>F)`[[1]])
    expect_true(grepl("median", lev_char$method))
})

test_that("Re-factor groups (useful to redo levels)", {
    notLow <- groups != "low"
    groups <- groups[notLow]
    values <- values[notLow]
    
    lev <- leveneTest(values, groups)
    car <- car::leveneTest(values, groups)
    
    expect_is(lev, "htest")
    expect_equal(lev$statistic[[1]], car$`F value`[[1]])
    expect_equal(lev$p.value, car$`Pr(>F)`[[1]])
    expect_true(grepl("median", lev$method))
})

test_that("Named values are calculated just as unnamed values", {
    names(values) <- paste0("name", 1:length(values))
    
    lev <- leveneTest(values, groups)
    car <- car::leveneTest(values, groups)
    
    expect_is(lev, "htest")
    expect_equal(lev$statistic[[1]], car$`F value`[[1]])
    expect_equal(lev$p.value, car$`Pr(>F)`[[1]])
    expect_true(grepl("median", lev$method))
})
