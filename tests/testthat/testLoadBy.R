context("Hierarchy check")

test_that("Gives FALSE if no one's responsible", {
    expect_false(loadBy(print, "someone"))
})

test_that("Gives boolean value depending on the responsible party", {
    d <- e <- f <- g <- h <- print
    attr(d, "loader") <- "charlie"
    attr(e, "loader") <- "snoopy"
    attr(f, "loader") <- "charlie"
    attr(g, "loader") <- "great pumpkin"
    attr(h, "loader") <- "woodstock"
    
    # Functions that should be loaded by "charlie"
    expect_true (loadBy("charlie", d))
    expect_false(loadBy("charlie", e))
    expect_true (loadBy("charlie", f))
    expect_false(loadBy("charlie", g))
    expect_false(loadBy("charlie", h))
    
    # Functions that should be loaded by "snoopy"
    expect_false(loadBy("snoopy", d))
    expect_true (loadBy("snoopy", e))
    expect_false(loadBy("snoopy", f))
    expect_false(loadBy("snoopy", g))
    expect_false(loadBy("snoopy", h))
    
    # Functions that should be loaded by "woodstock"
    expect_false(loadBy("woodstock", d))
    expect_false(loadBy("woodstock", e))
    expect_false(loadBy("woodstock", f))
    expect_false(loadBy("woodstock", g))
    expect_true (loadBy("woodstock", h))
})
