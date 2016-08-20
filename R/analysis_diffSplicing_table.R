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

#' Interface for differential analyses on all splicing events
#' 
#' @param id Character: identifier
#' 
#' @importFrom shinyjs disabled
#' 
#' @return HTML elements
diffSplicingTableUI <- function(id) {
    ns <- NS(id)
    
    sidebar <- sidebarPanel(
        selectizeInput(ns("groupsCol"), choices=NULL,
                       "Clinical groups on which to perform the analyses"), 
        uiOutput(ns("groupsInfo")), hr(),
        checkboxGroupInput(
            ns("statsChoices"),
            "Choose statistical analyses to perform:",
            # Basic stats is on and disabled by JavaScript
            c("Variance and median"="basicStats",
              "Wilcoxon signed rank test (1 group)"="wilcoxSignedRank",
              "Wilcoxon rank sum test (2 groups)"="wilcoxRankSum",
              "Kruskal-Wallis rank sum test (2 or more groups)"="kruskal", 
              "Levene's test (2 or more groups)"="levene",
              "Alternative splicing quantification density"="density"),
            selected=c("basicStats", "kruskal", "levene", "density",
                       "wilcoxSignedRank", "wilcoxRankSum")),
        # Disable checkbox of basic statistics
        tags$script('$("[value=basicStats]").attr("disabled", true);'),
        helpText("For each alternative splicing event, groups with one or less",
                 "non-missing values will be discarded."),
        disabled(div(id=ns("downloadStats"), class="btn-group",
                     tags$button(class="btn btn-default dropdown-toggle",
                                 type="button", "data-toggle"="dropdown",
                                 "aria-haspopup"="true",
                                 "aria-expanded"="false", 
                                 icon("download"), 
                                 "Download table", tags$span(class="caret")),
                     tags$ul(class="dropdown-menu", 
                             tags$li(downloadLink(ns("downloadAll"), 
                                                  "All data")),
                             tags$li(downloadLink(ns("downloadSubset"), 
                                                  "Filtered data"))))),
        processButton(ns("startAnalyses"), "Perform analyses"),
        uiOutput(ns("survivalOptions"))
    )
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebar, 
            mainPanel(
                dataTableOutput(ns("statsTable")),
                highchartOutput(ns("densitySparklines"), 0, 0)
            )
        )
    )
}

#' Create density sparklines for inclusion levels
#' @param data Character: HTML-formatted data series of interest
#' @param events Characater: event identifiers
#' 
#' @importFrom highcharter highchart hc_credits hc_tooltip hc_chart hc_title
#' hc_xAxis hc_yAxis hc_exporting hc_legend hc_plotOptions
#' @importFrom jsonlite toJSON
#' 
#' @return HTML element with sparkline data (character)
createDensitySparklines <- function(data, events) {
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
    
    json <- paste0(hc, ',"series":[', data, "]}")
    sparklines <- sprintf(
        paste('<sparkline onclick="showDiffSplicing(\'%s\')"',
              'style="cursor:pointer;" data-sparkline=\'%s\'/>'), 
        events, json)
    return(sparklines)
}

#' Interface for calculating optimal cut-off and p-value for survival curves
#' differences
#' @param ns Namespace function
optimSurvDiffUI <- function(ns) {
    tagList(
        hr(),
        h3("Survival analyses by splicing quantification cut-off"),
        helpText("For each splicing event, find the optimal splicing",
                 "quantification cut-off that most significantly separates",
                 "survival curves."),
        radioButtons(ns("censoring"), "Data censoring", selected="right",
                     inline=TRUE, choices=c(Left="left", Right="right",
                                            Interval="interval", 
                                            "Interval 2"="interval2")),
        selectizeInput(ns("timeStart"), choices = NULL, "Follow up time"),
        # If the chosen censoring contains the word 'interval', show this input
        conditionalPanel(paste0(
            "input[id='", ns("censoring"), "'].indexOf('interval') > -1"),
            selectizeInput(ns("timeStop"), choices=NULL, "Ending time")),
        helpText("In case there's no record for a patient, the",
                 "days to last follow up will be used instead."),
        selectizeInput(ns("event"), choices = NULL, 
                       "Event of interest"),
        radioButtons(
            ns("selected"), "Perform survival analysis in:",
            choices=c("Alternative splicing events shown in the screen"="shown",
                      "All alternative splicing events (slow process)"="all")),
        actionButton(ns("survival"), class="btn-primary",
                     "Perform survival analysis")
    )
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
        levene <- tryCatch(levene.test(vector[!nas], group[!nas]), error=return)
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
#' @param psi Data frame or matrix: Alternative splicing event quantification
#' @param analyses Character: analyses to perform
#' @param progress Function to track the progress
#' 
#' @importFrom plyr rbind.fill
#' @importFrom fastmatch fmatch
#' 
#' @details 
#' The following statistical analyses may be performed by including the 
#' respective string in the \code{analysis} argument:
#' \itemize{
#'      \item{Wilcoxon Rank Sum test - \code{wilcoxRankSum}}
#'      \item{Wilcoxon Signed Rank test - \code{wilcoxSignedRank}}
#'      \item{Kruskal test - \code{kruskal}}
#'      \item{Levene's test - \code{levene}}
#'      \item{Density plots - \code{density} (only usable through the visual 
#'      interface)}
#' }
#' 
#' @return Table of statistical analyses
#' @export
statsAnalyses <- function(psi, groups=NULL, analyses=c("wilcoxRankSum",
                                                       "wilcoxSignedRank",
                                                       "kruskal", "levene"),
                          progress=printPaste) {
    # cl <- parallel::makeCluster(getOption("cl.cores", getCores()))
    step <- 50 # Avoid updating progress for too few events
    progress("Performing statistical analysis", 
             divisions=5 + round(nrow(psi)/step))
    time <- Sys.time()
    
    if (is.null(groups)) {
        ids <- names(psi)
        groups <- parseSampleGroups(ids)
    }
    
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
    ll <- list()
    for (k in seq_along(uniq)) {
        df <- lapply(stats[match == k], data.frame)
        ll <- c(ll, list(do.call(rbind, df)))
    }
    
    # Bind all data frames together
    df <- do.call(rbind.fill, ll)
    rownames(df) <- unlist(lapply(ll, rownames))
    print(Sys.time() - time)
    
    # Remove columns of no interest
    df <- df[, !grepl("method|data.name", colnames(df))]
    
    # Calculate delta variance and delta median if there are only 2 groups
    deltaVar <- df[, grepl("Variance", colnames(df)), drop=FALSE]
    if (ncol(deltaVar) == 2) {
        progress("Calculating delta variance and median")
        deltaVar <- deltaVar[, 2] - deltaVar[, 1]
        deltaMed <- df[, grepl("Median", colnames(df))]
        deltaMed <- deltaMed[, 2] - deltaMed[, 1]
        df <- cbind(df, deltaVar, deltaMed)
    }
    
    if (any("density" == analyses)) {
        progress("Calculating the density of inclusion levels")
        df[, "Density"] <- createDensitySparklines(df[, "Density"],
                                                   rownames(df))
    }
    # parallel::stopCluster(cl)
    print(Sys.time() - time)
    return(df)
}

#' Optimal survival difference given an inclusion level cut-off for a specific
#' alternative splicing event
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
optimSurvDiff <- function(session, input, output) {
    # Interface of survival analyses
    output$survivalOptions <- renderUI({
        if (is.null(getDifferentialAnalyses()) || 
            is.null(getClinicalData()))
            return(NULL)
        
        optimSurvDiffUI(session$ns)
    })
    
    # Update clinical parameters
    observe({
        if (is.null(getDifferentialAnalyses()) || 
            is.null(getClinicalData()))
            return(NULL)
        
        updateClinicalParams(session)
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        if (is.null(getDifferentialAnalyses()) || 
            is.null(getClinicalData()) || is.null(input$censoring)) 
            return(NULL)
        
        label <- "Follow up time"
        if (grepl("interval", input$censoring, fixed=TRUE))
            label <- "Starting time"
        updateSelectizeInput(session, "timeStart", label=label)
    })
    
    #' Calculate optimal survival cut-off for the inclusion levels of a given
    #' alternative splicing event
    observeEvent(input$survival, {
        survTerms <- isolate({
            # Get tumour sample IDs (normal and control samples are not
            # interesting for survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- parseSampleGroups(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            # Group samples by the inclusion levels cut-off
            clinical <- getClinicalData()
            clinicalIDs <- nrow(clinical)
            groups <- rep(NA, clinicalIDs)
            
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
            psi       <- getInclusionLevels()
            stats     <- getDifferentialAnalyses()
            display   <- input$statsTable_rows_current
            selected  <- input$selected
        })
        
        if (selected == "shown") {
            if (!is.null(display)) {
                psi <- psi[display, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it's not possible to get events",
                           "shown in the table. To calculate survival analyses",
                           "calculate for all events.")
                closeProgress()
                return(NULL)
            }
        }
        startProgress("Performing survival analysis", nrow(psi))
        
        opt <- apply(psi, 1, function(vector) {
            v <- as.numeric(vector[toupper(names(tumour))])
            
            opt <- suppressWarnings(
                optim(0, testSurvivalCutoff, data=v, filter=tumour,
                      group=groups, clinical=clinical, censoring=censoring,
                      timeStart=timeStart, timeStop=timeStop, event=event,
                      # Method and parameters interval
                      method="Brent", lower=0, upper=1))
            
            updateProgress("Survival analysis", console=FALSE)
            return(c("Optimal survival PSI cut-off"=opt$par,
                     "Minimal survival p-value"=opt$value))
        })
        
        if (length(opt) == 0) {
            errorModal(session, "No survival analyses",
                       "Optimal PSI cut-off for the selected alternative",
                       "splicing events returned no survival analyses.")
        } else {
            # Remove NAs and add information to the statistical table
            opt <- opt[ , !is.na(opt[2, ])]
            df <- data.frame(t(opt))
            for (col in names(df)) stats[rownames(df), col] <- df[ , col]
            
            setDifferentialAnalyses(stats)
        }
        closeProgress()
        
        infoModal(session, "Survival columns added to table",
                  "The optimal survival cut-off and associated p-value for the",
                  "requested alternative splicing events were placed in the",
                  "last two columns of the table.")
    })
}

#' Server logic of the exploratory differential analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shinyjs toggleState
#' @importFrom DT replaceData dataTableProxy
diffSplicingTableServer <- function(input, output, session) {
    ns <- session$ns
    
    # Information on the data groups from TCGA
    output$groupsInfo <- renderUI({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        
        if (is.null(psi)) return(tagList(
            helpText(icon("exclamation-circle"), 
                     "No alternative splicing quantification loaded.",
                     tags$a(href="#", "Load or calculate it.",
                            onclick=loadRequiredData("Inclusion levels")))))
        
        # Separate samples by their type
        ids <- names(psi)
        type <- parseSampleGroups(ids)
        
        bullet <- "\u2022"
        groups <- NULL
        for (each in unique(type))
            groups <- tagList(groups, br(), bullet, each)
        
        return(tagList(
            helpText("The data contains the following sample types:", groups)))
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    
    performStatsAnalyses <- reactive({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        col <- input$groupsCol
        statsChoices <- input$statsChoices
        
        startProcessButton("startAnalyses")
        if (col == "Sample types") {
            # Separate samples by their groups
            ids <- names(psi)
            groups <- parseSampleGroups(ids)
        } else {
            # Get groups from column of interest
            clinical <- getClinicalData()
            col <- clinical[[col]]
            
            # Match groups from patients with respective samples
            matches <- getClinicalMatchFrom("Inclusion levels")
            groups <- rep(NA, ncol(psi))
            names(groups) <- colnames(psi)
            samples <- toupper(names(matches))
            groups[samples] <- as.character(col[matches])
            
            # Remove samples with no groups
            nasGroups <- !is.na(groups)
            psi       <- psi[nasGroups]
            groups    <- groups[nasGroups]
        }
        
        stats <- statsAnalyses(psi, groups, statsChoices, 
                               progress=updateProgress)
        
        stats <- cbind(stats, "Optimal survival PSI cut-off"=NA,
                       "Minimal survival p-value"=NA)
        
        setDifferentialAnalyses(stats)
        closeProgress()
        endProcessButton("startAnalyses")
    })
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            col <- input$groupsCol
            statsChoices <- input$statsChoices
            diffSplicing <- getDifferentialAnalyses()
        })
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
        } else if (is.null(col)) {
            errorModal(session, "Select groups",
                       "The groups on which to perform statistical analysis",
                       "cannot be empty.")
        } else if (!is.null(diffSplicing)) {
            warningModal(session, "Differential analyses already performed",
                         "Do you wish to replace the loaded analyses?",
                         footer=actionButton(ns("replace"), "Replace",
                                             class="btn-warning",
                                             "data-dismiss"="modal"))
        } else {
            performStatsAnalyses()
        }
    })
    
    observeEvent(input$replace, performStatsAnalyses())
    
    # proxy <- dataTableProxy("statsTable")
    # observe({
    #     # Do not re-render whole table if only the survival data is added
    #     if (!is.null(input$survival) && input$survival > 0) {
    #         replaceData(proxy, getDifferentialAnalyses())
    #     } else {
    output$statsTable <- renderDataTableSparklines({
        getDifferentialAnalyses()
    }, style="bootstrap", selection="none", filter='top', server=TRUE,
    extensions="Buttons", options=list(
        pageLength=10, rowCallback=JS("createDiffSplicingLinks"), dom='Bfrtip',
        buttons=I('colvis'),
        columnDefs=list(list(targets=1, searchable=FALSE))))
    # }
    # })
    
    # Download whole table
    output$downloadAll <- downloadHandler(
        filename=paste(getCategories(), "Differential splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalyses()
            densityCol <- NULL
            if (!is.null(stats)) densityCol <- match("Density", colnames(stats))
            
            write.table(stats[-densityCol], file, quote=FALSE, sep="\t",
                        row.names=FALSE)
        }
    )
    
    # Download filtered table
    output$downloadSubset <- downloadHandler(
        filename=paste(getCategories(), "Differential splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalyses()
            densityCol <- NULL
            if (!is.null(stats)) densityCol <- match("Density", colnames(stats))
            
            write.table(stats[input$statsTable_rows_all, -densityCol], 
                        file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )
    
    # Optimal survival difference given an inclusion level cut-off for a 
    # specific alternative splicing event
    optimSurvDiff(session, input, output)
    
    # Disable download button if statistical table is NULL
    observe({
        if (is.null(getDifferentialAnalyses()))
            disable("downloadStats")
        else
            enable("downloadStats")
    })
    
    # Update groups columns
    observe({
        clinical <- getClinicalData()
        psi <- getInclusionLevels()
        
        if (!is.null(clinical)) {
            updateSelectizeInput(
                session, "groupsCol", choices=list(
                    "Clinical groups for samples"=c("Sample types"="Sample types"),
                    "Clinical groups for patients"=names(clinical),
                    "d"=c("Start typing to search for clinical groups"="")))
        } else if (!is.null(psi)) {
            updateSelectizeInput(
                session, "groupsCol", choices=list(
                    "Clinical groups for samples"=c("Sample types"="Sample types"),
                    "d"=c("No clinical data loaded"="")))
        } else {
            updateSelectizeInput(session, "groupsCol", 
                                 choices=c("No clinical data loaded"=""))
        }
    })
    
    # Update groups used for differential splicing analysis
    observe( setDiffSplicingGroups(input$groupsCol) )
}

attr(diffSplicingTableUI, "loader") <- "diffSplicing"
attr(diffSplicingTableUI, "name") <- "All events (table)"
attr(diffSplicingTableUI, "selectEvent") <- FALSE
attr(diffSplicingTableServer, "loader") <- "diffSplicing"