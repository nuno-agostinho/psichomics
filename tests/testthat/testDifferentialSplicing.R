context("Differential splicing analysis")

# Calculate PSI
eventType <- c("SE", "MXE")

# Prepare annotation for SE and MXE
annot <- NULL
annot[[eventType[1]]] <- read.table(text = "1 + 32 35 37 38
                                    2 + 32 35 37 38
                                    3 + 32 35 37 38")
names(annot[[eventType[1]]]) <- c("Chromosome", "Strand",
                                  "C1.end", "A1.start", "A1.end", "C2.start")
annot[[eventType[2]]] <- read.table(text = "1 + 32 35 37 38 40 42
                                    2 + 32 35 37 38 40 42
                                    3 + 32 35 37 38 40 42")
names(annot[[eventType[2]]]) <- c("Chromosome", "Strand",
                                  "C1.end", "A1.start", "A1.end",
                                  "A2.start", "A2.end", "C2.start")

# Prepare junction quantification
junctionQuant <- read.table(text = "10 10 10 10 10 10
                            10 10 10 10 10 10
                            10 10 10 10 10 10
                            10 10 10 10 10 10
                            10 10 10 10 10 10
                            27 20 90 24 14 35
                            30 24 92 26 13 29
                            30 24 92 26 13 29
                            10 18 13 12 10 21
                            10 14 12 16 13 19
                            10 18 13 12 10 21
                            10 14 12 16 13 19
                            27 20 90 24 14 35
                            30 24 92 26 13 29
                            30 24 92 26 13 29")
colnames(junctionQuant) <- c(paste("Normal", 1:3), paste("Cancer", 1:3))
rownames(junctionQuant) <- c("chr1:32:+,chr1:35:+",
                             "chr1:37:+,chr1:42:+",
                             "chr1:37:+,chr1:38:+",
                             "chr1:32:+,chr1:38:+",
                             "chr1:40:+,chr1:42:+",
                             "chr2:37:+,chr2:38:+",
                             "chr2:32:+,chr2:35:+",
                             "chr2:37:+,chr2:42:+",
                             "chr2:32:+,chr2:38:+",
                             "chr2:40:+,chr2:42:+",
                             "chr3:32:+,chr3:35:+",
                             "chr3:37:+,chr3:42:+",
                             "chr3:32:+,chr3:38:+",
                             "chr3:40:+,chr3:42:+",
                             "chr3:37:+,chr3:38:+")
psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
group <- c(rep("Normal", 3), rep("Tumour", 3))

test_that("Perform all statistical analyses", {
    analyses <- c("wilcoxRankSum", "wilcoxSignedRank", "kruskal", "levene")
    stats <- statsAnalyses(psi, group, analyses)
    
    expect_is(stats, "data.frame")
    expect_named(stats)
    expect_true(all(stats$Samples.Normal == 3))
    expect_true(all(stats$Samples.Tumour == 3))
    expect_true(any(grepl("Wilcox", names(stats))))
    expect_true(any(grepl("Levene", names(stats))))
    expect_true(any(grepl("Kruskal", names(stats))))
    expect_true(any(grepl("Median", names(stats))))
    expect_true(any(grepl("Variance", names(stats))))
})

test_that("Perform single statistical analyses", {
    analyses <- c("wilcoxRankSum", "wilcoxSignedRank")
    stats <- statsAnalyses(psi, group, analyses)
    
    expect_is(stats, "data.frame")
    expect_named(stats)
    expect_true(all(stats$Samples.Normal == 3))
    expect_true(all(stats$Samples.Tumour == 3))
    expect_true(any(grepl("Wilcox", names(stats))))
    expect_false(any(grepl("Levene", names(stats))))
    expect_false(any(grepl("Kruskal", names(stats))))
    expect_true(any(grepl("Median", names(stats))))
    expect_true(any(grepl("Variance", names(stats))))
    
    analyses <- c("kruskal")
    stats <- statsAnalyses(psi, group, analyses)
    
    expect_is(stats, "data.frame")
    expect_named(stats)
    expect_true(all(stats$Samples.Normal == 3))
    expect_true(all(stats$Samples.Tumour == 3))
    expect_false(any(grepl("Wilcox", names(stats))))
    expect_false(any(grepl("Levene", names(stats))))
    expect_true(any(grepl("Kruskal", names(stats))))
    expect_true(any(grepl("Median", names(stats))))
    expect_true(any(grepl("Variance", names(stats))))
    
    analyses <- c("levene")
    stats <- statsAnalyses(psi, group, analyses)
    
    expect_is(stats, "data.frame")
    expect_named(stats)
    expect_true(all(stats$Samples.Normal == 3))
    expect_true(all(stats$Samples.Tumour == 3))
    expect_false(any(grepl("Wilcox", names(stats))))
    expect_true(any(grepl("Levene", names(stats))))
    expect_false(any(grepl("Kruskal", names(stats))))
    expect_true(any(grepl("Median", names(stats))))
    expect_true(any(grepl("Variance", names(stats))))
})

test_that("Plot distribution of splicing quantification per group", {
    eventPSI <- as.numeric(psi[6, ])
    plot <- plotDistribution(eventPSI, group)
    
    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equivalent(plot$x$hc_opts$xAxis[c("min", "max")], c(0, 1))
    expect_equal(plot$x$hc_opts$chart$zoomType, "x")
    # Plot two data series for each group (one of them is the rug plot)
    expect_length(plot$x$hc_opts$series, length(unique(group))*2)
})