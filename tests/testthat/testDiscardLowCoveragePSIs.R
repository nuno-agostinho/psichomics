context("Discard low coverage PSI values")

createPSIdataset <- function(samples=20, events=50) {
    set.seed(639245)
    size <- samples * events
    num  <- sample(seq(100), size, replace=TRUE)
    num[sample(size, size / 10)] <- NA

    psi <- matrix(num, ncol=samples) / 100
    colnames(psi) <- paste0("sample-", seq(ncol(psi)))
    rownames(psi) <- paste0("event-", seq(nrow(psi)))
    psi <- data.frame(psi)

    # Prepare event data
    quality <- c("N", "VLOW", "LOW", "OK", "SOK")
    quality <- sample(quality, size, replace=TRUE)
    eventData <- sprintf("SOK,%s,10=10=0,OK,S@20.00,0.00", quality)

    eventData <- data.frame(matrix(eventData, ncol=ncol(psi)))
    eventData <- cbind("vast-tools", eventData)
    colnames(eventData) <- c("source", paste0(colnames(psi), "-Q"))
    rownames(eventData) <- rownames(psi)
    class(eventData) <- c("eventData", class(eventData))
    attr(psi, "rowData") <- eventData

    psi <- preserveAttributes(psi)
    return(psi)
}

hasCvgValue <- function(eventData, cvg) {
    pattern <- sprintf(".*,%s,.*,.*,.*,.*", cvg)
    sapply(eventData, function(i) grepl(pattern, i))
}

checkDiscardCvgPSIvalues <- function(vals, samples=100, events=60) {
    psi    <- createPSIdataset(samples=samples, events=events)
    filter <- discardLowCoveragePSIvalues(psi, vasttoolsScoresToDiscard=vals)

    eventData <- attr(psi, "rowData")
    eventData <- eventData[ , endsWith(colnames(eventData), "-Q")]

    toNA <- matrix(FALSE, ncol=samples, nrow=events)
    for (val in vals) toNA <- toNA | hasCvgValue(eventData, val)

    psi[toNA] <- NA
    # Remove rows of missing values
    psi  <- psi[rowSums(!is.na(psi)) > 0, ]
    toNA <- toNA[rowSums(!toNA) > 0, ]

    noNAs <- nrow(toNA) == 1 && !any(toNA)
    if (nrow(filter) == 0) {
        expect_true(nrow(filter) == 0)
    } else if (noNAs) {
        expect_true(noNAs)
    } else {
        expect_true( all(is.na(unique(filter[toNA]))) )
    }
    expect_equivalent(psi, filter)
}

test_that("Discard low coverage VAST-TOOLS' PSI values", {
    skip_on_bioc()
    # Test one coverage value
    checkDiscardCvgPSIvalues("N")
    checkDiscardCvgPSIvalues("VLOW")
    checkDiscardCvgPSIvalues("LOW")
    checkDiscardCvgPSIvalues("OK")
    checkDiscardCvgPSIvalues("SOK")

    # Test two or more coverage values
    checkDiscardCvgPSIvalues(c("N", "VLOW"))
    checkDiscardCvgPSIvalues(c("OK", "SOK"))
    checkDiscardCvgPSIvalues(c("N", "VLOW", "LOW", "OK", "SOK"))

    # Test different input's sample size
    vals <- c("N", "VLOW", "OK", "SOK")
    checkDiscardCvgPSIvalues(vals)
    checkDiscardCvgPSIvalues(vals, samples=4, events=2)
    checkDiscardCvgPSIvalues(vals, samples=20, events=15)
})
