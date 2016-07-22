context("Functions to perform and plot principal component analysis (PCA)")

data <- USArrests
test_that("Perform PCA", {
    pca <- performPCA(data, center=FALSE, scale.=FALSE)
    expect_is(pca, "prcomp")
    expect_equal(nrow(pca$x), nrow(data))
    expect_equal(ncol(pca$x), ncol(data))
    
    # Internal PCA calculation used is "prcomp"
    pca2 <- prcomp(data, center=FALSE, scale.=FALSE)
    expect_equal(pca, pca2)
})

test_that("Center and scale the data", {
    center <- FALSE
    scale  <- FALSE
    pca <- performPCA(data, center=center, scale.=scale, naTolerance=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
    
    center <- TRUE
    scale  <- FALSE
    pca <- performPCA(data, center=center, scale.=scale, naTolerance=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
    
    center <- FALSE
    scale  <- TRUE
    pca <- performPCA(data, center=center, scale.=scale, naTolerance=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
    
    center <- TRUE
    scale  <- TRUE
    pca <- performPCA(data, center=center, scale.=scale, naTolerance=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
})

test_that("Tolerate NAs per row", {
    # Data has no NAs (therefore, naTolerance is irrelevant)
    pca <- performPCA(data, center=FALSE, scale.=FALSE, naTolerance=50)
    expect_is(pca, "prcomp")
    expect_equal(nrow(pca$x), nrow(data))
    expect_equal(ncol(pca$x), ncol(data))
    pca2 <- prcomp(data, center=FALSE, scale.=FALSE)
    expect_equal(pca, pca2)
    
    # Data is exclusively composed of NAs
    all.nas <- matrix(ncol=4, nrow=50)
    pca <- performPCA(all.nas, center=FALSE, scale.=FALSE, naTolerance=30)
    expect_null(pca)
    expect_error(prcomp(all.nas, center=center, scale.=scale))
    
    # Data is composed of roughly 50% NAs
    nas.50 <- data
    nas.50[replicate(4, sample(1:2, 50, replace=TRUE)) == 1] <- NA
    
    # Tolerate rows containing 75% of NAs
    pca <- performPCA(nas.50, center=FALSE, scale.=FALSE, naTolerance=75)
    expect_is(pca, "prcomp")
    expect_equal(rownames(nas.50)[rowSums(is.na(nas.50)) <= 3], rownames(pca$x))
    
    # Tolerate rows containing 50% of NAs
    pca <- performPCA(nas.50, center=FALSE, scale.=FALSE, naTolerance=50)
    expect_is(pca, "prcomp")
    expect_equal(rownames(nas.50)[rowSums(is.na(nas.50)) <= 2], rownames(pca$x))
    
    # Tolerate rows containing 25% of NAs
    pca <- performPCA(nas.50, center=FALSE, scale.=FALSE, naTolerance=25)
    expect_is(pca, "prcomp")
    expect_equal(rownames(nas.50)[rowSums(is.na(nas.50)) <= 1], rownames(pca$x))
})

test_that("Plot explained variance", {
    pca <- performPCA(data, center=FALSE, scale.=FALSE, naTolerance=0)
    hc <- plotVariance(pca)
    expect_is(hc, "highchart")
    expect_equal(hc$x$hc_opts$series[[1]]$data, pca$sdev ^ 2)
})

# test_that("Plot PCA", {
#     pca <- performPCA(data, center=FALSE, scale.=FALSE)
#     imp <- summary(pca)$importance[2, ]
#     perc <- as.numeric(imp)
#     names(perc) <- names(imp)
#     
#     plotPCA(pca, perc, xAxis="PC1", yAxis="PC2", selected=rownames(pca$x))
# })