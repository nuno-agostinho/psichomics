context("Correlation between gene expression and PSI values")

# Calculate PSI
eventType <- c("SE", "MXE")
annot <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")
psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
rownames(psi) <- gsub("(?:.(?!_))+$", "_", rownames(psi), perl=TRUE)
rownames(psi) <- paste0(rownames(psi), c("A1K", "BDK3", "BDA1KL",
                                         "A1K", "MHN/A1K/BDK3", "XHR/MHNOR"))
geneExpr <- data.frame(c(20, 51, 32, 50, 60, 90),
                       c(10, 40, 30, 20, 50, 60),
                       c(96, 45, 96, 65, 24, 12),
                       c(63, 56, 72, 89, 83, 27),
                       c(30, 59, 74, 12, 25, 32),
                       c(36, 77, 41, 54, 12, 58))
rownames(geneExpr) <- c("A1K|423", "BDK3|567", "BDA1KL|754", "MHN|2120", 
                        "MHNOR|2134", "XHR|12442")
group <- c(rep("Normal", 3), rep("Cancer", 3))
colnames(geneExpr) <- paste(group, 1:3)

test_that("Plot the correlation between GE and AS quantification", {
    expect_warning(corr <- correlateGEandAS(geneExpr, psi, "A1K"))
    expect_equal(names(corr), rownames(psi)[c(1, 4, 5)])
    expect_length(plotCorrelation(corr), 3)
    
    corr <- correlateGEandAS(geneExpr, psi, "BDK3")
    expect_equal(names(corr), rownames(psi)[c(2, 5)])
    expect_length(plotCorrelation(corr), 2)
    
    corr <- correlateGEandAS(geneExpr, psi, "XHR")
    expect_equal(names(corr), rownames(psi)[6])
    expect_length(plotCorrelation(corr), 1)
})

# Test with NAs within the PSI values
psi[4, 2] <- psi[2, 6] <- psi[1, 3] <- psi[5, 6] <- psi[3, 3] <- psi[5, 4] <- NA

test_that("Correlate gene expression and AS quantification with NAs", {
    expect_warning(corr <- correlateGEandAS(geneExpr, psi, "A1K"))
    expect_equal(names(corr), rownames(psi)[c(1, 4, 5)])
    expect_length(plotCorrelation(corr), 3)
    
    corr <- correlateGEandAS(geneExpr, psi, "BDK3")
    expect_equal(names(corr), rownames(psi)[c(2, 5)])
    expect_length(plotCorrelation(corr), 2)
    
    corr <- correlateGEandAS(geneExpr, psi, "MHN")
    expect_equal(names(corr), rownames(psi)[5])
    expect_length(plotCorrelation(corr), 1)
})

test_that("Correctly match genes based on TCGA-styled gene expression", {
    # Test using a mixed of partial and full matching genes
    gene <- c("A1K|423", "BDK3|567", "BDA1KL", "MHN|2120", "MHNOR", "XHR|12442")
    subset <- subsetGeneExpressionFromMatchingGenes(geneExpr, gene)
    expect_equal(rownames(subset),
                 c("A1K|423", "BDK3|567", "BDA1KL|754", "MHN|2120", 
                   "MHNOR|2134", "XHR|12442"))
    
    # Test with only TCGA-styled genes (full match)
    gene <- c("A1K|423", "BDK3|567", "MHN|2120", "XHR|12442")
    subset <- subsetGeneExpressionFromMatchingGenes(geneExpr, gene)
    expect_equal(rownames(subset), 
                 c("A1K|423", "BDK3|567", "MHN|2120", "XHR|12442"))
    
    # Test with only non-TCGA-styled genes (partial match)
    gene <- c("BDA1KL", "MHNOR")
    subset <- subsetGeneExpressionFromMatchingGenes(geneExpr, gene)
    expect_equal(rownames(subset), c("BDA1KL|754", "MHNOR|2134"))
    
    # Test with non-matching genes
    gene   <- c("SOMESPAM", "A1K|423", "BDK3", "MORESPAM", "MHN", "SPAMSPAM")
    subset <- subsetGeneExpressionFromMatchingGenes(geneExpr, gene)
    expect_equal(rownames(subset), c("A1K|423", "BDK3|567", "MHN|2120"))
    
    # Test with only non-matching genes
    gene <- c("ALLSPAM", "SPAMSAMPLE", "SOMEMORESPAM")
    expect_error(subsetGeneExpressionFromMatchingGenes(geneExpr, gene))
})