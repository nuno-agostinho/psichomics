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

test_that("rowMeans calculates the mean per row of a matrix", {
    # Passing a matrix
    mat <- replicate(10, rnorm(20))
    precisionError <- 10e-16
    test <- apply(mat, 1, mean) - rowMeans(mat) < precisionError
    expect_true(all(test))
    
    # Passing a single vector
    mat <- mat[1, ]
    expect_equal(rowMeans(mat), mean(mat))
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
    
    match <- getPatientFromSample(samples, clinical)
    expect_is(match, "character")
    expect_equivalent(match, patients)
    expect_equal(names(match), samples)
})

test_that("Retrieve samples from patient identifiers", {
    patients <- c("GTEX-ABC", "GTEX-DEF", "GTEX-GHI", "GTEX-JKL", "GTEX-MNO")
    samples <- paste0(patients, "-sample")
    clinical <- data.frame(samples=samples)
    rownames(clinical) <- patients
    
    ref <- c(1, 4)
    match <- getMatchingSamples(patients[ref], samples, clinical)
    expect_is(match, "character")
    expect_equivalent(match[], as.character(clinical$samples[ref]))
    
    # Retrieve samples when previously matched
    patients <- getPatientFromSample(samples, clinical)
    match <- getMatchingSamples(patients[ref], samples, clinical, 
                                match=patients)
    expect_is(match, "character")
    expect_equivalent(match[], as.character(clinical$samples[ref]))
})

test_that("Parse alternative splicing event from identifiers", {
    events <- c(
        "A3SS_15_+_63353138_63353912_63353397_TPM1",
        "A3SS_11_-_61118463_61117115_61117894_CYB561A3",
        "A5SS_21_+_48055675_48056459_48056808_PRMT2", 
        "A5SS_1_-_1274742_1274667_1274033_DVL1",
        "AFE_9_+_131902430_131901928_131904724_PPP2R4",
        "AFE_5_-_134686513_134688636_134681747_H2AFY",
        "ALE_12_+_56554104_56554410_56555171_MYL6",
        "ALE_8_-_38314874_38287466_38285953_FGFR1",
        "SE_9_+_6486925_6492303_6492401_6493826_UHRF2",
        "SE_19_-_5218431_5216778_5216731_5215606_PTPRS",
        "MXE_15_+_63335142_63335905_63336030_63336226_63336351_63349184_TPM1",
        "MXE_17_-_74090495_74087316_74087224_74086478_74086410_74085401_EXOC7")
    
    parsed <- parseSplicingEvent(events, coords=TRUE)
    expect_is(parsed, "data.frame")
    expect_equal(unique(parsed$type),
                 c("A3SS", "A5SS", "AFE", "ALE", "SE", "MXE"))
    expect_equal(as.numeric(parsed$chrom), 
                 c(15, 11, 21, 1, 9, 5, 12, 8, 9, 19, 15, 17))
    expect_equal(unlist(parsed$gene),
                 c("TPM1", "CYB561A3", "PRMT2", "DVL1", "PPP2R4", "H2AFY",
                   "MYL6", "FGFR1", "UHRF2", "PTPRS", "TPM1", "EXOC7"))
    expect_equal(tail(colnames(parsed), 4), 
                 c("constitutive1", "alternative1", "alternative2",
                   "constitutive2"))
})

