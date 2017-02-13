context("Test miscellanious functions")

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

test_that("Text suggestions create a runnable JS script", {
    words <- c("gene", "transcript", "protein")
    suggest <- textSuggestions("id", words)
    expect_is(suggest, "html")
    
    # Words are in script
    scriptWords <- paste0('["', paste(words, collapse = '", "'), "\"]")
    expect_true(grepl(scriptWords, suggest, fixed=TRUE))
    
    # The library textcomplete is used
    expect_true(grepl(".textcomplete(", suggest, fixed=TRUE))
})

test_that("Create a button with a loading indicator", {
    id <- "buttonId"
    label <- "Click me!"
    button <- processButton(id, label)
    expect_equal(button[[2]]$id, id)
    expect_equal(button[[2]]$type, "button")
    expect_equal(button[[3]][[1]][[2]][[3]][[2]], label)
    
    icon <- button[[3]][[1]][[2]][[3]][[1]]
    expect_equal(icon[[1]], "i")
    expect_equal(icon[[2]][[2]], "fa fa-spinner fa-spin")
    expect_equal(icon[[2]][[3]], "shinyjs-hide")
})

test_that("Retrieve patients from sample identifiers", {
    patients <- c("GTEX-ABC", "GTEX-DEF", "GTEX-GHI", "GTEX-JKL", "GTEX-MNO")
    samples <- paste0(patients, "-sample")
    clinical <- data.frame(samples=samples)
    rownames(clinical) <- patients
    
    patients <- getPatientFromSample(samples, clinical)
    expect_is(patients, "integer")
    expect_equivalent(patients, 1:5)
    expect_equal(names(patients), samples)
    
    # # Do not remove non-matching identifiers (by default)
    # clinical[4:5, "samples"] <- NA
    # patients <- getPatientFromSample(samples, clinical)
    # expect_equivalent(patients, c(1:3, NA, NA))
    # expect_equal(names(patients), samples)
    
    # # Remove non-matching identifiers
    # patients <- getPatientFromSample(samples, clinical)
    # expect_equivalent(patients, 1:3)
    # expect_equal(names(patients), samples[1:3])
})

test_that("Retrieve samples from patient identifiers", {
    patients <- c("GTEX-ABC", "GTEX-DEF", "GTEX-GHI", "GTEX-JKL", "GTEX-MNO")
    samples <- paste0(patients, "-sample")
    clinical <- data.frame(samples=samples)
    rownames(clinical) <- patients
    
    ref <- c(1, 4)
    match <- getMatchingSamples(ref, samples, clinical)
    expect_is(match, "character")
    expect_equivalent(match[], as.character(clinical$samples[ref]))
    
    # Retrieve samples when previously matched
    patients <- getPatientFromSample(samples, clinical)
    match <- getMatchingSamples(ref, samples, clinical, match=patients)
    expect_is(match, "character")
    expect_equivalent(match[], as.character(clinical$samples[ref]))
})