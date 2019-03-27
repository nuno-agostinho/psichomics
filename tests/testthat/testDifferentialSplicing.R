context("Differential splicing analysis")

# Calculate PSI
eventType <- c("SE", "MXE")
annot     <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")
psi       <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
group     <- c(rep("Normal", 3), rep("Tumour", 3))
group_ls  <- split(colnames(junctionQuant), group)

test_that("Perform all statistical analyses", {
    analyses <- c("wilcoxRankSum", "wilcoxSignedRank", "kruskal", "levene")
    stats <- diffAnalyses(psi, group, analyses)
    
    expect_is(stats, "data.frame")
    expect_named(stats)
    expect_true(all(stats$Samples.Normal == 3))
    expect_true(all(stats$Samples.Tumour == 3))
    expect_true(any(grepl("Wilcox", names(stats))))
    expect_true(any(grepl("Levene", names(stats))))
    expect_true(any(grepl("Kruskal", names(stats))))
    expect_true(any(grepl("Median", names(stats))))
    expect_true(any(grepl("Variance", names(stats))))

    # Perform with a list of groups
    analyses <- c("wilcoxRankSum", "wilcoxSignedRank", "kruskal", "levene")
    stats2 <- diffAnalyses(psi, group_ls, analyses)
    expect_identical(stats, stats2)
})

test_that("Perform all statistical analyses for a single event", {
    analyses <- c("wilcoxRankSum", "wilcoxSignedRank", "kruskal", "levene")
    stats <- diffAnalyses(psi[1, ], group, analyses)
    
    expect_is(stats, "data.frame")
    expect_named(stats)
    expect_true(all(stats$Samples.Normal == 3))
    expect_true(all(stats$Samples.Tumour == 3))
    expect_true(any(grepl("Wilcox", names(stats))))
    expect_true(any(grepl("Levene", names(stats))))
    expect_true(any(grepl("Kruskal", names(stats))))
    expect_true(any(grepl("Median", names(stats))))
    expect_true(any(grepl("Variance", names(stats))))
    
    stats2 <- diffAnalyses(psi[1, ], group_ls, analyses)
    expect_identical(stats, stats2)
})

test_that("Perform single statistical tests", {
    analyses <- c("wilcoxRankSum", "wilcoxSignedRank")
    stats <- diffAnalyses(psi, group, analyses)
    
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
    stats <- diffAnalyses(psi, group, analyses)
    
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
    stats <- diffAnalyses(psi, group, analyses)
    
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
    eventPSI <- as.numeric(psi[2, ])
    plot <- plotDistribution(eventPSI, group)
    
    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equivalent(plot$x$hc_opts$xAxis[c("min", "max")], c(0, 1))
    expect_equal(plot$x$hc_opts$chart$zoomType, "x")
    # Plot two data series for each group (one of them is the rug plot)
    expect_length(plot$x$hc_opts$series, length(unique(group))*2)
})

test_that("Label groups based on a cutoff", {
    data <- c(1, 0, 0, 1, 0.5, 1)
    cutoff <- 0.5
    
    # Greater or equal than a cutoff (default)
    group <- labelBasedOnCutoff(data, cutoff)
    group <- gsub("&gt;", ">", group)
    group <- gsub("&lt;", "<", group)
    expect_is(group, "character")
    expect_equal(group, paste(c(">=", "<", "<", ">=", ">=", ">="), cutoff))
    
    # Greater than a cutoff
    group <- labelBasedOnCutoff(data, cutoff, gte=FALSE)
    group <- gsub("&gt;", ">", group)
    group <- gsub("&lt;", "<", group)
    expect_equal(group, paste(c(">", "<=", "<=", ">", "<=", ">"), cutoff))
    
    # Add text to label before
    label <- "Proportion"
    group <- labelBasedOnCutoff(data, cutoff, label=label)
    group <- gsub("&gt;", ">", group)
    group <- gsub("&lt;", "<", group)
    expect_equal(group,
                 paste(label, c(">=", "<", "<", ">=", ">=", ">="), cutoff))
})
