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
    pca <- performPCA(data, center=center, scale.=scale, missingValues=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
    
    center <- TRUE
    scale  <- FALSE
    pca <- performPCA(data, center=center, scale.=scale, missingValues=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
    
    center <- FALSE
    scale  <- TRUE
    pca <- performPCA(data, center=center, scale.=scale, missingValues=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
    
    center <- TRUE
    scale  <- TRUE
    pca <- performPCA(data, center=center, scale.=scale, missingValues=0)
    pca2 <- prcomp(data, center=center, scale.=scale)
    expect_equal(pca, pca2)
})

test_that("Tolerate NAs per columns", {
    # Data has no NAs (therefore, missingValues is irrelevant)
    pca <- performPCA(data, center=FALSE, scale.=FALSE, missingValues=50)
    expect_is(pca, "prcomp")
    expect_equal(nrow(pca$x), nrow(data))
    expect_equal(ncol(pca$x), ncol(data))
    pca2 <- prcomp(data, center=FALSE, scale.=FALSE)
    expect_equal(pca, pca2)
    
    # Data is exclusively composed of NAs
    all.nas <- matrix(ncol=4, nrow=50)
    expect_warning(
        pca <- performPCA(all.nas, center=FALSE, scale.=FALSE, missingValues=30))
    expect_null(pca)
    expect_error(prcomp(all.nas, center=FALSE, scale.=FALSE))
    
    # Fill with missing values (column 1 = 100% NAs, 2 = 50%, 3 = 34%, 4 = 26%)
    nas <- data
    nas[[1]][seq(1, length(nas[[1]]), 1)] <- NA
    nas[[2]][seq(1, length(nas[[2]]), 2)] <- NA
    nas[[3]][seq(1, length(nas[[3]]), 3)] <- NA
    nas[[4]][seq(1, length(nas[[4]]), 4)] <- NA
    
    # Tolerate columns containing 50% of NAs
    pca <- performPCA(nas, center=FALSE, scale.=FALSE, missingValues=25)
    expect_equal(colnames(nas)[2:4], rownames(pca$rotation))
    
    # Tolerate columns containing 49% of NAs
    pca <- performPCA(nas, center=FALSE, scale.=FALSE, missingValues=24)
    expect_equal(colnames(nas)[3:4], rownames(pca$rotation))
    
    # Tolerate columns containing 26% of NAs
    pca <- performPCA(nas, center=FALSE, scale.=FALSE, missingValues=13)
    expect_equal(colnames(nas)[4], rownames(pca$rotation))
    
    # Tolerate columns containing 25% of NAs
    expect_warning(
        pca <- performPCA(nas, center=FALSE, scale.=FALSE, missingValues=12))
    expect_null(pca)
})

test_that("Plot explained variance", {
    pca <- performPCA(data, center=FALSE, scale.=FALSE, missingValues=0)
    hc <- plotVariance(pca)
    expect_is(hc, "highchart")
    eigenvalue <- vapply(hc$x$hc_opts$series[[1]]$data, "[[", "eigenvalue", 
                         FUN.VALUE = numeric(1))
    expect_equal(eigenvalue, pca$sdev ^ 2)
})

pca <- performPCA(data, center=FALSE, scale.=FALSE)
groups <- lapply(1:4, `*`, c(2, 5))
names(groups) <- paste("Group", 1:4)
groups <- lapply(groups, function(i) rownames(pca$x)[i])

test_that("Plot all PCA individuals", {
    hc <- plotPCA(pca, pcX="PC1", pcY="PC2")
    expect_is(hc, "highchart")
    
    opts <- hc$x$hc_opts
    expect_null(sapply(opts$series, "[[", "name")[[1]])
})

test_that("Plot PCA individuals and colour all groups", {
    hc <- plotPCA(pca, pcX="PC1", pcY="PC2", groups)
    expect_is(hc, "highchart")
    
    opts <- hc$x$hc_opts
    expect_equal(sapply(opts$series, "[[", "name"), names(groups))
})

test_that("Plot PCA individuals and colour two groups", {
    hc <- plotPCA(pca, pcX="PC1", pcY="PC2", groups[2:3])
    expect_is(hc, "highchart")
    
    opts <- hc$x$hc_opts
    expect_equal(sapply(opts$series, "[[", "name"), names(groups)[2:3])
})

test_that("Plot PCA loadings", {
    hc <- plotPCA(pca, pcX="PC1", pcY="PC2", loadings=TRUE)
    expect_is(hc, "highchart")
    
    opts <- hc$x$hc_opts
    expect_is(opts$series[[2]], "list")
    
    # Colour two groups of individuals
    hc <- plotPCA(pca, pcX="PC1", pcY="PC2", groups[2:3], loadings=TRUE)
    expect_is(hc, "highchart")
    
    opts <- hc$x$hc_opts
    namz <- sapply(opts$series, "[[", "name")
    expect_equal(unlist(namz), c(names(groups)[2:3], "Loadings"))

    # Plot different principal components
    hc <- plotPCA(pca, pcX="PC3", pcY="PC4", loadings=TRUE)
    expect_is(hc, "highchart")
    expect_equal(hc$x$hc_opts$series[[2]]$name, "Loadings")
    expect_equal(hc$x$hc_opts$series[[2]]$type, "bubble")
})
