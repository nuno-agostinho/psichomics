#' Modified from lawstat::levene.test
#' @inheritParams lawstat::levene.test
#' @importFrom stats anova lm
#' @importFrom stats kruskal.test weighted.mean
levene.test <- function (y, group, location = c("median", "mean", "trim.mean"), 
                         trim.alpha = 0.25, bootstrap = FALSE, num.bootstrap = 1000, 
                         kruskal.test = FALSE, correction.method = c("none", "correction.factor", 
                                                                     "zero.removal", "zero.correction")) 
{
    if (length(y) != length(group)) {
        stop("the length of the data (y) does not match the length of the group")
    }
    location <- match.arg(location)
    correction.method <- match.arg(correction.method)
    DNAME = deparse(substitute(y))
    y <- y[!is.na(y)]
    group <- group[!is.na(y)]
    if ((location == "trim.mean") & (trim.alpha == 1)) {
        stop("trim.alpha value of 0 to 0.5 should be provided for the trim.mean location")
    }
    reorder <- order(group)
    group <- group[reorder]
    y <- y[reorder]
    gr <- group
    group <- as.factor(group)
    if (location == "mean") {
        means <- tapply(y, group, mean)
        METHOD <- "classical Levene's test based on the absolute deviations from the mean"
    }
    else if (location == "median") {
        means <- tapply(y, group, median)
        METHOD = "modified robust Brown-Forsythe Levene-type test based on the absolute deviations from the median"
    }
    else {
        location = "trim.mean"
        trimmed.mean <- function(y) mean(y, trim = trim.alpha)
        means <- tapply(y, group, trimmed.mean)
        METHOD <- "modified robust Levene-type test based on the absolute deviations from the trimmed mean"
    }
    n <- tapply(y, group, length)
    resp.mean <- abs(y - means[group])
    ngroup <- n[group]
    if (location != "median" && correction.method != "correction.factor") {
        METHOD <- paste(METHOD, "(", correction.method, "not applied because the location is not set to median", 
                        ")")
        correction.method <- "none"
    }
    if (correction.method == "correction.factor") {
        METHOD <- paste(METHOD, "with correction factor")
        correction <- sqrt(ngroup/(ngroup - 1))
        resp.mean <- correction * resp.mean
    }
    if (correction.method == "zero.removal" || correction.method == 
        "zero.correction") {
        if (correction.method == "zero.removal") {
            METHOD <- paste(METHOD, "with Hines-Hines structural zero removal method")
        }
        if (correction.method == "zero.correction") {
            METHOD <- paste(METHOD, "with modified structural zero removal method and correction factor")
        }
        resp.mean <- y - means[group]
        k <- length(n)
        temp <- double()
        endpos <- double()
        startpos <- double()
        for (i in 1:k) {
            group.size <- n[i]
            j <- i - 1
            if (i == 1) 
                start <- 1
            else start <- sum(n[1:j]) + 1
            startpos <- c(startpos, start)
            end <- sum(n[1:i])
            endpos <- c(endpos, end)
            sub.resp.mean <- resp.mean[start:end]
            sub.resp.mean <- sub.resp.mean[order(sub.resp.mean)]
            if (group.size%%2 == 1) {
                mid <- (group.size + 1)/2
                temp2 <- sub.resp.mean[-mid]
                if (correction.method == "zero.correction") {
                    ntemp <- length(temp2) + 1
                    correction <- sqrt((ntemp - 1)/ntemp)
                    temp2 <- correction * temp2
                }
            }
            if (group.size%%2 == 0) {
                mid <- group.size/2
                if (correction.method == "zero.removal") {
                    denom <- sqrt(2)
                }
                else {
                    denom <- 1
                }
                replace1 <- (sub.resp.mean[mid + 1] - sub.resp.mean[mid])/denom
                temp2 <- sub.resp.mean[c(-mid, -mid - 1)]
                temp2 <- c(temp2, replace1)
                if (correction.method == "zero.correction") {
                    ntemp <- length(temp2) + 1
                    correction <- sqrt((ntemp - 1)/ntemp)
                    temp2 <- correction * temp2
                }
            }
            temp <- c(temp, temp2)
        }
        resp.mean <- abs(temp)
        zero.removal.group <- group[-endpos]
    }
    else {
        correction.method = "none"
    }
    if (correction.method == "zero.removal" || correction.method == 
        "zero.correction") {
        d <- zero.removal.group
    }
    else {
        d <- group
    }
    if (kruskal.test == FALSE) {
        nv <- anova(lm(resp.mean ~ d))
        statistic <- nv[1, 4]
        p.value <- nv[1, 5]
    }
    else {
        METHOD <- paste("rank-based (Kruskal-Wallis)", METHOD)
        ktest <- kruskal.test(resp.mean, d)
        statistic <- ktest$statistic
        p.value = ktest$p.value
    }
    non.bootstrap.p.value <- p.value
    if (bootstrap == TRUE) {
        METHOD <- paste("bootstrap", METHOD)
        R <- 0
        N <- length(y)
        frac.trim.alpha <- 0.2
        b.trimmed.mean <- function(y) {
            nn <- length(y)
            wt <- rep(0, nn)
            y2 <- y[order(y)]
            lower <- ceiling(nn * frac.trim.alpha) + 1
            upper <- floor(nn * (1 - frac.trim.alpha))
            if (lower > upper) 
                stop("frac.trim.alpha value is too large")
            m <- upper - lower + 1
            frac <- (nn * (1 - 2 * frac.trim.alpha) - m)/2
            wt[lower - 1] <- frac
            wt[upper + 1] <- frac
            wt[lower:upper] <- 1
            return(weighted.mean(y2, wt))
        }
        b.trim.means <- tapply(y, group, b.trimmed.mean)
        rm <- y - b.trim.means[group]
        for (j in 1:num.bootstrap) {
            sam <- sample(rm, replace = TRUE)
            boot.sample <- sam
            if (min(n) < 10) {
                U <- runif(1) - 0.5
                means <- tapply(y, group, mean)
                v <- sqrt(sum((y - means[group])^2)/N)
                boot.sample <- ((12/13)^(0.5)) * (sam + v * U)
            }
            if (location == "mean") {
                boot.means <- tapply(boot.sample, group, mean)
            }
            else if (location == "median") {
                boot.means <- tapply(boot.sample, group, median)
            }
            else {
                location = "trim.mean"
                trimmed.mean.2 <- function(boot.sample) mean(boot.sample, 
                                                             trim = trim.alpha)
                boot.means <- tapply(boot.sample, group, trimmed.mean.2)
            }
            resp.boot.mean <- abs(boot.sample - boot.means[group])
            if (correction.method == "correction.factor") {
                correction <- sqrt(ngroup/(ngroup - 1))
                resp.mean <- correction * resp.boot.mean
            }
            if (correction.method == "zero.removal" || correction.method == 
                "zero.correction") {
                resp.mean <- boot.sample - boot.means[group]
                k <- length(n)
                temp <- double()
                endpos <- double()
                startpos <- double()
                for (i in 1:k) {
                    group.size <- n[i]
                    j <- i - 1
                    if (i == 1) 
                        start <- 1
                    else start <- sum(n[1:j]) + 1
                    startpos <- c(startpos, start)
                    end <- sum(n[1:i])
                    endpos <- c(endpos, end)
                    sub.resp.mean <- resp.mean[start:end]
                    sub.resp.mean <- sub.resp.mean[order(sub.resp.mean)]
                    if (group.size%%2 == 1) {
                        mid <- (group.size + 1)/2
                        temp2 <- sub.resp.mean[-mid]
                        if (correction.method == "zero.correction") {
                            ntemp <- length(temp2) + 1
                            correction <- sqrt((ntemp - 1)/ntemp)
                            temp2 <- correction * temp2
                        }
                    }
                    if (group.size%%2 == 0) {
                        mid <- group.size/2
                        if (correction.method == "zero.removal") {
                            denom <- sqrt(2)
                        }
                        else {
                            denom <- 1
                        }
                        replace1 <- (sub.resp.mean[mid + 1] - sub.resp.mean[mid])/denom
                        temp2 <- sub.resp.mean[c(-mid, -mid - 1)]
                        temp2 <- c(temp2, replace1)
                        if (correction.method == "zero.correction") {
                            ntemp <- length(temp2) + 1
                            correction <- sqrt((ntemp - 1)/ntemp)
                            temp2 <- correction * temp2
                        }
                    }
                    temp <- c(temp, temp2)
                }
                resp.boot.mean <- abs(temp)
                zero.removal.group <- group[-endpos]
            }
            if (correction.method == "zero.removal" || correction.method == 
                "zero.correction") {
                d <- zero.removal.group
            }
            else {
                d <- group
            }
            if (kruskal.test == FALSE) {
                statistic2 = anova(lm(resp.boot.mean ~ d))[1, 
                                                           4]
            }
            else {
                bktest <- kruskal.test(resp.boot.mean, d)
                statistic2 <- bktest$statistic
            }
            if (statistic2 > statistic) 
                R <- R + 1
        }
        p.value <- R/num.bootstrap
    }
    STATISTIC = statistic
    names(STATISTIC) = "Test Statistic"
    structure(list(statistic = STATISTIC, p.value = p.value, 
                   method = METHOD, data.name = DNAME, non.bootstrap.p.value = non.bootstrap.p.value), 
              class = "htest")
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
#' @param events Characater: event identifiers
#' @param delim Character: left and right delimeters in groups that should be
#' removed
#' 
#' @importFrom highcharter highchart hc_credits hc_tooltip hc_chart hc_title
#' hc_xAxis hc_yAxis hc_exporting hc_legend hc_plotOptions
#' @importFrom jsonlite toJSON
#' 
#' @return HTML element with sparkline data (character)
createDensitySparklines <- function(data, events, delim=NULL) {
    hc <- highchart() %>%
        hc_tooltip(
            hideDelay=0, shared=TRUE, valueDecimals=getPrecision(),
            headerFormat="<small>Inclusion levels: {point.x}</small><br/>",
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
#'      \item{Wilcoxon Rank Sum test - \code{wilcoxRankSum}}
#'      \item{Wilcoxon Signed Rank test - \code{wilcoxSignedRank}}
#'      \item{Kruskal test - \code{kruskal}}
#'      \item{Levene's test - \code{levene}}
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
#' @importFrom stats kruskal.test median wilcox.test var density
#' 
#' @return A row from a data frame with the results
singleStatsAnalyses <- function(vector, group, threshold=1, step=100,
                                analyses=c("wilcoxRankSum", "wilcoxSignedRank",
                                           "kruskal", "levene")) {
    series  <- split(vector, group)
    samples <- vapply(series, function(i) sum(!is.na(i)), integer(1))
    valid   <- names(series)[samples >= threshold]
    inGroup <- group %in% valid
    group   <- group[inGroup]
    vector  <- vector[inGroup]
    len     <- length(valid)
    
    # Wilcoxon tests
    wilcox <- NULL
    if (any("wilcoxRankSum" == analyses) && len == 2) {
        # Wilcoxon rank sum test (Mann Whitney U test)
        typeOne <- group == valid[1]
        wilcox  <- suppressWarnings(wilcox.test(vector[typeOne],
                                                vector[!typeOne]))
    } else if (any("wilcoxSignedRank" == analyses) && len == 1) {
        # Wilcoxon signed rank test
        wilcox <- suppressWarnings(wilcox.test(vector))
    }
    
    # Kruskal-Wallis test
    kruskal <- NULL
    if (any("kruskal" == analyses) && len >= 2) {
        kruskal <- tryCatch(kruskal.test(vector, group), error=return)
        if (any("error" == class(kruskal))) kruskal <- NULL
    }
    
    # Levene's test
    levene <- NULL
    if (any("levene" == analyses) && len >= 2) {
        nas <- is.na(vector)
        levene <- suppressWarnings(
            tryCatch(levene.test(vector[!nas], group[!nas]), error=return))
        if (any("error" == class(levene))) levene <- NULL
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
    
    vector <- c(Density=sparkline, Samples=samples, Wilcox=wilcox, 
                Kruskal=kruskal, Levene=levene, Variance=var, Median=med)
    vector <- vector[!vapply(vector, is.null, logical(1))] # Remove NULL
    return(vector)
}

#' Perform selected statistical analyses on multiple splicing events
#' 
#' @param psi Data frame or matrix: alternative splicing event quantification
#' @param analyses Character: analyses to perform (see Details)
#' @param pvalueAdjust Character: method used to adjust p-values (see Details)
#' @param progress Function to track the progress
#' 
#' @importFrom plyr rbind.fill
#' @importFrom fastmatch fmatch
#' 
#' @details 
#' The following statistical analyses may be performed by including the 
#' respective string in the \code{analysis} argument:
#' \itemize{
#'      \item{\code{wilcoxRankSum}: Wilcoxon Rank Sum test}
#'      \item{\code{wilcoxSignedRank}: Wilcoxon Signed Rank test}
#'      \item{\code{kruskal}: Kruskal test}
#'      \item{\code{levene}: Levene's test}
#'      \item{\code{density}: Density plots (only usable through the visual 
#'      interface)}
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
statsAnalyses <- function(psi, groups=NULL, analyses=c("wilcoxRankSum",
                                                       "wilcoxSignedRank",
                                                       "kruskal", "levene"),
                          pvalueAdjust="BH", progress=printPaste) {
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
        return(singleStatsAnalyses(...))
    }, factor(groups), threshold=1, step=step, analyses=analyses)
    print(Sys.time() - time)
    
    # Check the column names of the different columns
    ns <- lapply(stats, names)
    uniq <- unique(ns)
    match <- fmatch(ns, uniq)
    
    # Convert matrix with the same name to data frames (way faster)
    progress("Preparing data")
    time <- Sys.time()
    ll <- list()
    for (k in seq_along(uniq)) {
        df <- lapply(stats[match == k], data.frame)
        ll <- c(ll, list(do.call(rbind, df)))
    }
    
    # Bind all data frames together
    df <- do.call(rbind.fill, ll)
    rownames(df) <- unlist(lapply(ll, rownames))
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
    
    if (any("density" == analyses)) {
        progress("Calculating the density of inclusion levels")
        time <- Sys.time()
        df[, "Density"] <- createDensitySparklines(
            df[, "Density"], rownames(df), delim=c(parenthesisOpen, 
                                                   parenthesisClose))
        print(Sys.time() - time)
    }
    
    if (any(pvalueAdjust == c("BH", "BY", "bonferroni", "holm", "hochberg",
                              "hommel"))) {
        progress("Adjusting p-values", detail=pvalueAdjust)
        time <- Sys.time()
        
        cols   <- grep("p.value", colnames(df), fixed=TRUE)
        pvalue <- df[cols]
        adjust <- apply(pvalue, 2, p.adjust, pvalueAdjust)
        names  <- paste0(colnames(pvalue), " ", parenthesisOpen, 
                         pvalueAdjust, " adjusted", parenthesisClose)
        
        if (!is.matrix(adjust))
            adjust <- as.matrix(adjust)
        colnames(adjust) <- names
        
        # Place the adjusted p-values columns next to the respective p-values
        len <- ncol(df)
        order <- seq(len)
        for (i in seq_along(cols)) 
            order <- append(order, len+i, after=which(order == cols[i]))
        df <- cbind(df, adjust)[order]
        print(Sys.time() - time)
    }
    
    # Properly set column names
    col <- colnames(df)
    col <- gsub(parenthesisOpen,  "(", col, fixed=TRUE)
    col <- gsub(parenthesisClose, ")", col, fixed=TRUE)
    col <- gsub(".", " ", col, fixed=TRUE)
    col <- gsub("p value", "p-value", col, fixed=TRUE)
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
diffSplicingServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("diffSplicing",
                                 priority=c("diffSplicingTableServer",
                                            "diffSplicingEventServer"))
}

attr(diffSplicingUI, "loader") <- "analysis"
attr(diffSplicingUI, "name") <- "Differential splicing analysis"
attr(diffSplicingServer, "loader") <- "analysis"