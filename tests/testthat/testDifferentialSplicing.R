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

context("Visualise sample distribution via density plots")
eventPSI <- as.numeric(psi[2, ])

test_that("Plot distribution of multiple values per group", {
    plot <- suppressWarnings(plotDistribution(eventPSI, group))

    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equal(plot$x$hc_opts$chart$zoomType, "x")
    expect_identical(sort(unique(sapply(plot$x$hc_opts$series, "[[", "name"))),
                     sort(unique(group)))

    # Plot two data series for each group: density (area) + rug (scatter) plot
    nGroups <- length(unique(group))
    expect_length(plot$x$hc_opts$series, nGroups * 2)
    expect_identical(sapply(plot$x$hc_opts$series, "[[", "type"),
                     rep(c("area", "scatter"), each=nGroups))

    # If values are within 0 and 1, assume data to be PSI values
    expect_equivalent(plot$x$hc_opts$xAxis[c("min", "max")], c(0, 1))
    expect_equal(plot$x$hc_opts$xAxis$title$text, "Distribution of PSI values")

    # If values are outside the range [0, 1], assume data to be gene expression
    value <- c(2, 5, 2, 4, 5, 6)
    plot <- plotDistribution(value, group)
    expect_null(plot$x$hc_opts$xAxis$min)
    expect_null(plot$x$hc_opts$xAxis$max)
    expect_equal(plot$x$hc_opts$xAxis$title$text,
                 "Distribution of gene expression")

    # No group assigned, i.e. all data values correspond to the same group
    plot  <- plotDistribution(value)
    expect_length(plot$x$hc_opts$series, 2) # Two series: density + rug plot
    expect_identical(plot$x$hc_opts$series[[1]]$name, "All samples")
    expect_identical(plot$x$hc_opts$series[[2]]$name, "All samples")
})

test_that("Plot distribution based on different types of groups", {
    # If values are named, groups can also be a list of names
    names(eventPSI) <- paste("Sample", seq(eventPSI))
    group2          <- split(names(eventPSI), group)
    p1 <- plotDistribution(eventPSI, group)
    p2 <- plotDistribution(eventPSI, group2)

    # Y values are randomly assigned: use same Y values to compare the two
    scatter <- which(sapply(p1$x$hc_opts$series, "[[", "type") == "scatter")
    for (i in scatter) {
        data <- p1$x$hc_opts$series[[i]]$data
        for (j in seq(data)) {
            p2$x$hc_opts$series[[i]]$data[[j]]$y <- data[[j]]$y
        }
    }
    expect_identical(p1, p2)
})

test_that("Plot distribution of a single value", {
    value <- 2.3
    plot  <- plotDistribution(value)
    series <- plot$x$hc_opts$series
    expect_length(series, 2)

    # Density (area) plot is invisible and disabled
    expect_equal(series[[1]]$type, "area")
    expect_false(series[[1]]$visible)
    expect_equivalent(as.character(series[[1]]$events$legendItemClick),
                      "function(e) { e.preventDefault() }")

    expect_equal(series[[2]]$type, "scatter")
    expect_equal(series[[2]]$data[[1]]$x, value)
})

test_that("Plot distribution can set custom title", {
    plot <- plotDistribution(eventPSI, group)
    expect_null(plot$x$hc_opts$title$text)

    title <- "Custom title"
    plot <- plotDistribution(eventPSI, group, title=title)
    expect_equal(plot$x$hc_opts$title$text, title)
})

test_that("Plot distribution with labels", {
    areRugLabelsEnabled <- function(plot) {
        type <- sapply(plot$x$hc_opts$series, "[[", "type") == "scatter"
        areLabelsEnabled <- sapply(plot$x$hc_opts$series[type],
                                   function(x) x[["dataLabels"]][["enabled"]])
        return(all(areLabelsEnabled))
    }

    areLabelsEquivalent <- function(plot, name, comp) {
        series <- unlist(plot$x$hc_opts$series)
        labels <- series[names(series) == paste0("data.", name)]
        expect_equivalent(labels, comp)
    }

    # If names are unset and rugLabels = FALSE: do not show rug + tooltip labels
    names(eventPSI) <- NULL
    plot <- plotDistribution(eventPSI, group, rugLabels=FALSE)
    expect_false(areRugLabelsEnabled(plot))
    areLabelsEquivalent(plot, "tooltipLabel", character(0))
    areLabelsEquivalent(plot, "tooltipLabel", character(0))
    isTooltipCorrectlyFormatted <- grepl("^\\{point\\.tooltipLabel\\}",
                                         plot$x$hc_opts$tooltip$pointFormat)
    expect_true(isTooltipCorrectlyFormatted)

    # If names are unset and rugLabels = TRUE: show rug labels only
    names(eventPSI) <- NULL
    plot <- plotDistribution(eventPSI, group, rugLabels=TRUE)
    expect_true(areRugLabelsEnabled(plot))
    areLabelsEquivalent(plot, "rugLabel", as.character(round(eventPSI, 2)))
    areLabelsEquivalent(plot, "tooltipLabel", character(0))

    # If names are set and rugLabels = FALSE: show tooltip labels only
    labels <- paste("Sample", seq(eventPSI))
    names(eventPSI) <- labels
    plot <- plotDistribution(eventPSI, group, rugLabels=FALSE)
    expect_false(areRugLabelsEnabled(plot))
    areLabelsEquivalent(plot, "rugLabel", character(0))
    areLabelsEquivalent(plot, "tooltipLabel", labels)

    # If names are set and rugLabels = TRUE: show both rug and tooltip labels
    labels <- paste("Sample", seq(eventPSI))
    names(eventPSI) <- labels
    plot <- plotDistribution(eventPSI, group, rugLabels=TRUE)
    expect_true(areRugLabelsEnabled(plot))
    areLabelsEquivalent(plot, "rugLabel", labels)
    areLabelsEquivalent(plot, "tooltipLabel", labels)
})

context("Label groups based on a cutoff")

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
