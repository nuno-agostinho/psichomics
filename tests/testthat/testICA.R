context("Functions to perform and plot independent component analysis (ICA)")

data <- USArrests
test_that("Perform ICA", {
    ica <- performICA(data, center=FALSE, scale.=FALSE)
    expect_is(ica, "list")
    expect_equal(nrow(ica$S), nrow(data))
    expect_equal(ncol(ica$S), ncol(data))
})

test_that("Tolerate NAs per columns", {
    # Data is exclusively composed of NAs
    all.nas <- matrix(ncol=4, nrow=50)
    expect_warning(
        ica <- performICA(all.nas, center=FALSE, scale.=FALSE, 
                          missingValues=30))
    expect_null(ica)
    
    # Fill with missing values (column 1 = 100% NAs, 2 = 50%, 3 = 34%, 4 = 26%)
    nas <- data
    nas[[1]][seq(1, length(nas[[1]]), 1)] <- NA
    nas[[2]][seq(1, length(nas[[2]]), 2)] <- NA
    nas[[3]][seq(1, length(nas[[3]]), 3)] <- NA
    nas[[4]][seq(1, length(nas[[4]]), 4)] <- NA
    
    # Tolerate columns containing 50% of NAs
    ica <- performICA(nas, center=FALSE, scale.=FALSE, missingValues=25)
    expect_equal(colnames(nas)[2:4], colnames(ica$X))
    
    # Tolerate columns containing 49% of NAs
    ica <- performICA(nas, center=FALSE, scale.=FALSE, missingValues=24)
    expect_equal(colnames(nas)[3:4], colnames(ica$X))
    
    # "Tolerate" columns containing 26% of NAs
    ica <- performICA(nas, center=FALSE, scale.=FALSE, missingValues=13)
    expect_is(ica, "error")
    
    # Tolerate columns containing 25% of NAs
    expect_warning(
        ica <- performICA(nas, center=FALSE, scale.=FALSE, missingValues=12))
    expect_null(ica)
})

ica <- performICA(data, center=FALSE, scale.=FALSE)
groups <- lapply(1:4, `*`, c(2, 5))
names(groups) <- paste("Group", 1:4)
groups <- lapply(groups, function(i) rownames(ica$S)[i])

test_that("Plot all ICA individuals", {
    hc <- plotICA(ica)
    expect_is(hc, "pairsD3")
    expect_true(all(ica$S == hc$x$data))
})

test_that("Plot ICA individuals and colour all groups", {
    hc <- plotICA(ica, 1:2, groups)
    expect_is(hc, "pairsD3")
    subset <- ica$S[unlist(groups), 1:2]
    expect_true(all(subset == hc$x$data))
})

test_that("Plot ICA individuals and colour two groups", {
    hc <- plotICA(ica, 1:2, groups[2:3])
    expect_is(hc, "pairsD3")
    subset <- ica$S[unlist(groups[2:3]), 1:2]
    expect_true(all(subset == hc$x$data))
})
