#' Levene's test
#' 
#' Performs a Levene's test to assess the equality of variances
#' 
#' @inheritParams stats::kruskal.test
#' @param centers Function used to calculate how much values spread 
#' (\code{median} by default; another common function used is \code{mean})
#' 
#' @importFrom stats complete.cases anova median lm
#' 
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic with a name describing it.}
#' \item{p.value}{the p-value for the test.}
#' \item{method}{the type of test applied.}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @examples 
#' 
#' vals <- sample(30, replace=TRUE)
#' group <- lapply(list("A", "B", "C"), rep, 10)
#' group <- unlist(group)
#' psichomics:::leveneTest(vals, group)
#' 
#' ## Using Levene's test based on the mean
#' psichomics:::leveneTest(vals, group, mean)
leveneTest <- function(x, g, centers=median) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    
    # Remove missing values
    noNAs <- complete.cases(x, g)
    x <- x[noNAs]
    g <- g[noNAs]
    
    # Convert groups to factors
    g <- factor(g)
    
    res <- vapply(split(x, g, drop=TRUE), centers, numeric(1))
    spread <- abs(x - res[g])
    
    # Analysis of variance (ANOVA)
    var <- anova(lm(spread ~ g))
    statistic <- var$`F value`[1]
    pval <- var$`Pr(>F)`[1]
    
    centers <- deparse(substitute(centers))
    rval <- list(statistic=c("W"=statistic), p.value=pval, data.name=dname,
                 method=paste0("Levene's test (using the ", centers, ")"))
    class(rval) <- "htest"
    return(rval)
}

#' User interface for the differential splicing analyses
#' 
#' @param id Character: identifier
#' @param tab Function to process HTML elements
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS div icon fluidRow column selectizeInput conditionalPanel
#' 
#' @return HTML element as character
diffSplicingUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "diffSplicing", 
                             priority=c("diffSplicingTableUI",
                                        "diffSplicingEventUI"))
    
    ui <- lapply(uiList, function(ui) tabPanel(attr(ui, "name"), ui) )
    do.call(tabsetPanel, c(list(type="pills"), ui))
}

#' Create density sparklines for inclusion levels
#' @param data Character: HTML-formatted data series of interest
#' @param events Character: event identifiers
#' @param delim Character: left and right delimeters in groups that should be
#' removed
#' 
#' @importFrom highcharter highchart hc_credits hc_tooltip hc_chart hc_title
#' hc_xAxis hc_yAxis hc_exporting hc_legend hc_plotOptions
#' @importFrom jsonlite toJSON
#' @importFrom shiny tags
#' 
#' @return HTML element with sparkline data (character)
createDensitySparklines <- function(data, events, delim=NULL) {
    hc <- highchart() %>%
        hc_tooltip(
            hideDelay=0, shared=TRUE, valueDecimals=getPrecision(),
            headerFormat=paste0(tags$small("Inclusion levels: {point.x:.2f}"),
                                tags$br()),
            pointFormat=paste(span(style="color:{point.color}", "\u25CF "),
                              tags$b("{series.name}"), br())) %>%
        hc_chart(width=120, height=20, backgroundColor="", type="areaspline", 
                 margin=c(2, 0, 2, 0), style=list(overflow='visible')) %>%
        hc_title(text="") %>%
        hc_credits(enabled=FALSE) %>%
        hc_xAxis(min=0, max=1, visible=FALSE) %>%
        hc_yAxis(endOnTick=FALSE, startOnTick=FALSE, visible=FALSE) %>%
        hc_exporting(enabled=FALSE) %>%
        hc_legend(enabled=FALSE) %>%
        hc_plotOptions(series=list(cursor="non", animation=FALSE, lineWidth=1,
                                   marker=list(radius=1), fillOpacity=0.25))
    hc <- as.character(toJSON(hc$x$hc_opts, auto_unbox=TRUE))
    hc <- substr(hc, 1, nchar(hc)-1)
    
    # Remove artificial delimeters
    if (!is.null(delim)) {
        data <- gsub(delim[[1]], "", data)
        data <- gsub(delim[[2]], "", data)
    }
    
    json <- paste0(hc, ',"series":[', data, "]}")
    sparklines <- sprintf(
        paste('<sparkline onclick="showDiffSplicing(\'%s\')"',
              'style="cursor:pointer;" data-sparkline=\'%s\'/>'), 
        events, json)
    return(sparklines)
}

#' Perform statistical analysis on a given splicing event
#' 
#' Perform statistical analyses on a given vector containing elements from
#' different groups
#' 
#' @details 
#' The following statistical analyses may be performed by including the 
#' respective string in the \code{analysis} argument:
#' \itemize{
#'      \item{\code{ttest} - Unpaired t-test (2 groups)}
#'      \item{\code{wilcoxRankSum} - Wilcoxon Rank Sum test (2 groups)}
#'      \item{\code{kruskal} - Kruskal test (2 or more groups)}
#'      \item{\code{levene} - Levene's test (2 or more groups)}
#'      \item{\code{fligner} - Fligner-Killeen test (2 or more groups)}
#' }
#' 
#' @param vector Numeric
#' @param group Character: group of each element in the vector
#' @param threshold Integer: minimum number of data points to perform analysis
#' in a group (default is 1)
#' @param analyses Character: analyses to perform (see "Details")
#' @param step Numeric: number of events before the progress bar is updated
#' (a bigger number allows for a faster execution)
#' 
#' @importFrom stats kruskal.test median wilcox.test t.test var density
#' @importFrom methods is
#' 
#' @return A row from a data frame with the results
singleDiffAnalyses <- function(vector, group, threshold=1, step=100,
                               analyses=c("wilcoxRankSum", "ttest", "kruskal",
                                          "levene", "fligner")) {
    series  <- split(vector, group)
    samples <- vapply(series, function(i) sum(!is.na(i)), integer(1))
    valid   <- names(series)[samples >= threshold]
    inGroup <- group %in% valid
    group   <- group[inGroup]
    vector  <- vector[inGroup]
    len     <- length(valid)
    
    # Unpaired t-test (2 groups)
    ttest <- NULL
    if (any("ttest" == analyses) && len == 2 && all(samples > 1)) {
        typeOne <- group == valid[1]
        ttest <- tryCatch(t.test(vector[typeOne], vector[!typeOne]), 
                          error=return)
        if (is(ttest, "error")) ttest <- NULL
    }
    
    # Wilcoxon test (2 groups)
    wilcox <- NULL
    if (any("wilcoxRankSum" == analyses) && len == 2) {
        # Wilcoxon rank sum test (Mann Whitney U test)
        typeOne <- group == valid[1]
        wilcox  <- suppressWarnings(wilcox.test(vector[typeOne],
                                                vector[!typeOne]))
    # } else if (any("wilcoxSignedRank" == analyses) && len == 1) {
    #     # Wilcoxon signed rank test
    #     wilcox <- suppressWarnings(wilcox.test(vector))
    }
    
    # Kruskal-Wallis test (2 or more groups)
    kruskal <- NULL
    if (any("kruskal" == analyses) && len >= 2) {
        kruskal <- tryCatch(kruskal.test(vector, group), error=return)
        if (any("error" == class(kruskal))) kruskal <- NULL
    }
    
    # Levene's test (2 or more groups)
    levene <- NULL
    if (any("levene" == analyses) && len >= 2) {
        levene <- suppressWarnings(
            tryCatch(leveneTest(vector, group), error=return))
        if (any("error" == class(levene))) levene <- NULL
    }
    
    # Fligner-Killeen test (2 or more groups)
    fligner <- NULL
    if (any("fligner" == analyses) && len >= 2) {
        fligner <- suppressWarnings(
            tryCatch(fligner.test(vector, group), error=return))
        if (any("error" == class(fligner))) fligner <- NULL
    }
    
    # Density sparklines
    sparkline <- NULL
    if (any("density" == analyses)) {
        data <- NULL
        validSeries <- series[valid]
        for (group in validSeries) {
            # Calculate the density of inclusion levels for each sample type 
            # with a greatly reduced number of points for faster execution
            den <- density(group, n=10, bw=0.01, na.rm=TRUE)
            data <- c(data, paste(sprintf('{"x":%s,"y":%s}', den$x, den$y),
                                  collapse=","))
        }
        sparkline <- paste(sprintf('{"name":"%s", "data":[%s]}', 
                                   names(validSeries), data), collapse=",")
    }
    
    # Variance and median
    med <- lapply(series, median, na.rm=TRUE) # Median
    var <- lapply(series, var, na.rm=TRUE) # Variance
    
    vector <- c("PSI.distribution"=sparkline, Samples=samples, "T-test"=ttest, 
                Wilcoxon=wilcox, Kruskal=kruskal, Levene=levene, 
                "Fligner-Killeen"=fligner, Variance=var, Median=med)
    vector <- vector[!vapply(vector, is.null, logical(1))] # Remove NULL
    return(vector)
}

#' Perform selected statistical analyses on multiple splicing events
#' 
#' @param psi Data frame or matrix: alternative splicing event quantification
#' @param groups Character: group of each sample from the alternative splicing 
#' event quantification (if NULL, sample types are used instead, e.g. normal, 
#' tumour and metastasis)
#' @param analyses Character: analyses to perform (see Details)
#' @param pvalueAdjust Character: method used to adjust p-values (see Details)
#' @param progress Function to track the progress
#' 
#' @importFrom plyr rbind.fill
#' @importFrom fastmatch fmatch
#' @importFrom stats p.adjust
#' 
#' @details 
#' The following statistical analyses may be performed by including the 
#' respective string in the \code{analysis} argument:
#' \itemize{
#'      \item{\code{ttest} - Unpaired t-test (2 groups)}
#'      \item{\code{wilcoxRankSum} - Wilcoxon Rank Sum test (2 groups)}
#'      \item{\code{kruskal} - Kruskal test (2 or more groups)}
#'      \item{\code{levene} - Levene's test (2 or more groups)}
#'      \item{\code{fligner} - Fligner-Killeen test (2 or more groups)}
#'      \item{\code{density} - Sample distribution per group (only usable 
#'      through the visual interface)}
#' }
#' 
#' The following methods for p-value adjustment are supported by using the 
#' respective string in the \code{pvalueAdjust} argument:
#' \itemize{
#'      \item{\code{none}: do not adjust p-values}
#'      \item{\code{BH}: Benjamini-Hochberg's method (false discovery rate)}
#'      \item{\code{BY}: Benjamini-Yekutieli's method (false discovery rate)}
#'      \item{\code{bonferroni}: Bonferroni correction (family-wise error rate)}
#'      \item{\code{holm}: Holm's method (family-wise error rate)}
#'      \item{\code{hochberg}: Hochberg's method (family-wise error rate)}
#'      \item{\code{hommel}: Hommel's method (family-wise error rate)}
#' }
#' 
#' @return Table of statistical analyses
#' @export
#' @examples 
#' # Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
#' eventType <- c("SE", "MXE")
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' 
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#' group <- c(rep("Normal", 3), rep("Tumour", 3))
#' diffAnalyses(psi, group)
diffAnalyses <- function(psi, groups=NULL, 
                         analyses=c("wilcoxRankSum", "ttest", "kruskal",
                                    "levene", "fligner"),
                         pvalueAdjust="BH", progress=echoProgress) {
    # cl <- parallel::makeCluster(getOption("cl.cores", getCores()))
    step <- 50 # Avoid updating progress too frequently
    progress("Performing statistical analysis", 
             divisions=5 + round(nrow(psi)/step))
    time <- Sys.time()
    
    if (is.null(groups)) {
        ids <- names(psi)
        groups <- parseSampleGroups(ids)
    }
    
    # Add artificial delimeters (required to identify group names later on)
    parenthesisOpen  <- ".delim1."
    parenthesisClose <- ".delim2."
    groups <- paste0(parenthesisOpen, groups, parenthesisClose)
    
    count <- 0
    stats <- apply(psi, 1, function(...) {
        count <<- count + 1
        if (count %% step == 0)
            progress("Performing statistical analysis", console=FALSE)
        return(singleDiffAnalyses(...))
    }, factor(groups), threshold=1, step=step, analyses=analyses)
    print(Sys.time() - time)
    
    # Check the column names of the different columns
    ns <- lapply(stats, names)
    uniq <- unique(ns)
    match <- fmatch(ns, uniq)
    
    progress("Preparing data")
    time <- Sys.time()

    # Convert list of lists to data frame
    ll <- lapply(stats, function(i) lapply(i, unname))
    ll <- lapply(ll, unlist)
    ldf <- lapply(seq_along(uniq), function(k) {
        elems <- match == k
        df2 <- data.frame(t(as.data.frame(ll[elems])))
        rownames(df2) <- names(stats)[elems]
        return(df2)
    })
    df <- do.call(rbind.fill, ldf)
    rownames(df) <- unlist(lapply(ldf, rownames))
    
    # Convert numeric columns to numeric
    num <- suppressWarnings(apply(df, 2, as.numeric))
    if (!is.matrix(num)) {
        num <- t(as.matrix(num))
        rownames(num) <- rownames(df)
    }
    numericCols <- colSums(is.na(num)) != nrow(num)
    df[ , numericCols] <- num[ , numericCols]
    
    # Convert integer columns to integer
    if (any(numericCols)) {
        int <- apply(df[ , numericCols, drop=FALSE], 2, function(i) 
            all(is.whole(i), na.rm=TRUE))
        intCols <- numericCols
        intCols[numericCols] <- int
        if (any(intCols))
            df[ , intCols] <- apply(df[ , intCols, drop=FALSE], 2, as.integer)
    }
    print(Sys.time() - time)
    
    # Calculate delta variance and delta median if there are only 2 groups
    deltaVar <- df[, grepl("Variance", colnames(df)), drop=FALSE]
    if (ncol(deltaVar) == 2) {
        progress("Calculating delta variance and median")
        time <- Sys.time()
        deltaVar <- deltaVar[, 2] - deltaVar[, 1]
        deltaMed <- df[, grepl("Median", colnames(df))]
        deltaMed <- deltaMed[, 2] - deltaMed[, 1]
        df <- cbind(df, "Delta variance"=deltaVar, "Delta median"=deltaMed)
        print(Sys.time() - time)
    }
    
    if (any(pvalueAdjust == c("BH", "BY", "bonferroni", "holm", "hochberg",
                              "hommel"))) {
        progress("Adjusting p-values", detail=pvalueAdjust)
        
        cols   <- grep("p.value", colnames(df), fixed=TRUE)
        if (length(cols > 0)) {
            time <- Sys.time()
            pvalue <- df[cols]
            adjust <- apply(pvalue, 2, p.adjust, pvalueAdjust)
            names  <- paste0(colnames(pvalue), " ", parenthesisOpen, 
                             pvalueAdjust, " adjusted", parenthesisClose)
            
            if (!is.matrix(adjust))
                adjust <- t(as.matrix(adjust))
            colnames(adjust) <- names
            
            # Place the adjusted p-values next to the respective p-values
            len <- ncol(df)
            order <- seq(len)
            for (i in seq_along(cols)) 
                order <- append(order, len+i, after=which(order == cols[i]))
            df <- cbind(df, adjust)[order]
            print(Sys.time() - time)
        }
    }
    
    # Add splicing event information
    progress("Include splicing event information")
    info <- suppressWarnings(parseSplicingEvent(rownames(df)))
    
    # Prepare presentation of multigenes
    multigene <- lapply(info$gene, length) > 1
    infoGene <- info$gene
    infoGene[multigene] <- lapply(infoGene[multigene], paste, collapse="/")
    
    df <- cbind("Event type"=info$type, "Chromosome"=info$chrom,
                "Strand"=info$strand, "Gene"=unlist(infoGene), df)
    
    if (any("density" == analyses)) {
        progress("Calculating the density of inclusion levels")
        time <- Sys.time()
        df[, "PSI.distribution"] <- createDensitySparklines(
            df[, "PSI.distribution"], rownames(df), delim=c(parenthesisOpen, 
                                                   parenthesisClose))
        print(Sys.time() - time)
    }
    
    # Properly set column names
    col <- colnames(df)
    col <- gsub(parenthesisOpen,  "(", col, fixed=TRUE)
    col <- gsub(parenthesisClose, ")", col, fixed=TRUE)
    col <- gsub(".", " ", col, fixed=TRUE)
    col <- gsub("p value", "p-value", col, fixed=TRUE)
    col <- gsub("T test", "T-test", col, fixed=TRUE)
    col <- gsub("Fligner Killeen", "Fligner-Killeen", col, fixed=TRUE)
    colnames(df) <- col
    
    # parallel::stopCluster(cl)
    return(df)
}

#' Server logic for the differential splicing analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny observe observeEvent updateSelectizeInput
#' @importFrom shinyjs hide show
#' @return NULL (this function is used to modify the Shiny session's state)
diffSplicingServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("diffSplicing",
                                 priority=c("diffSplicingTableServer",
                                            "diffSplicingEventServer"))
}

attr(diffSplicingUI, "loader") <- "analysis"
attr(diffSplicingUI, "name") <- "Differential splicing analysis"
attr(diffSplicingServer, "loader") <- "analysis"