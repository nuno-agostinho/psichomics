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

#' Interface for exploratory differential analyses
#' @param id Character: identifier
#' @return HTML elements
diffAnalysisTableUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebarPanel(
                uiOutput(ns("groupsInfo")),
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
                helpText("If groups have too many missing values for a",
                         "particular alternative splicing event, the group",
                         "itself will be discarded from the analyses."),
                actionButton(ns("startAnalyses"), class="btn-primary", 
                             "Perform analyses"),
                uiOutput(ns("survivalOptions"))
            ), mainPanel(
                uiOutput(ns("showColumns")),
                dataTableOutput(ns("statsTable")),
                highchartOutput(ns("densitySparklines"), 0, 0)
            )
        )
    )
}

#' Perform statistical analysis on a vector with elements from different groups
#' 
#' @param vector Numeric
#' @param group Character: group of each element in the vector
#' @param threshold Integer: minimum number of data points to perform analysis
#' in a group (default is 1)
#' @param analyses Character: name of the analyses to perform (all by default)
#' @param step Numeric: number of events before the progress bar is updated
#' (a bigger number allows for a faster execution)
#' 
#' @importFrom stats kruskal.test median wilcox.test var
#' 
#' @return A data frame row with the results
statsAnalyses <- function(vector, group, threshold=1, step=100,
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

#' Calculate density sparklines for inclusion levels
#' @param data Character: HTML-formatted data series of interest
#' 
#' @importFrom highcharter highchart hc_credits hc_tooltip hc_chart hc_title
#' hc_xAxis hc_yAxis hc_exporting hc_legend hc_plotOptions
#' @importFrom jsonlite toJSON
#' 
#' @return HTML element with sparkline data (character)
calculateDensitySparklines <- function(data) {
    hc <- highchart() %>%
        hc_tooltip(
            hideDelay=0, shared=TRUE,
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
    sparklines <- sprintf("<sparkline data-sparkline='%s'/>", json)
    return(sparklines)
}

#' Interface for calculating optimal cut-off and p-value for survival curves
#' differences
#' @param ns Namespace function
optimSurvDiffUI <- function(ns) {
    tagList(
        hr(),
        h3("Survival analyses"),
        radioButtons(ns("censoring"), "Data censoring",
                     selected="right",
                     inline=TRUE, choices=c(
                         Left="left",
                         Right="right",
                         Interval="interval",
                         "Interval 2" = "interval2")),
        selectizeInput(ns("timeStart"), choices = NULL, 
                       "Follow up time"),
        # If the chosen censoring contains the word 'interval',
        # show this input
        conditionalPanel(
            paste0("input[id='", ns("censoring"),
                   "'].indexOf('interval') > -1"),
            selectizeInput(ns("timeStop"), choices=NULL, 
                           "Ending time")),
        helpText("In case there's no record for a patient, the",
                 "days to last follow up will be used instead."),
        selectizeInput(ns("event"), choices = NULL, 
                       "Event of interest"),
        radioButtons(
            ns("selected"), "Perform survival analysis in:",
            choices=c("Alternative splicing events currently in the table"="shown",
                      "All alternative splicing events (slow process)"="all")),
        actionButton(ns("survival"), class="btn-primary",
                     "Perform survival analysis")
    )
}

#' Perform all statistical analyses
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param psi Data frame or matrix: Alternative splicing event quantification
#' @param statsChoices Character: analyses to perform
#' 
#' @importFrom plyr rbind.fill
#' @importFrom fastmatch fmatch
performStatsAnalyses <- function(session, input, psi, statsChoices) {
    # Separate samples by their type
    ids <- names(psi)
    type <- getSampleTypes(ids)
    
    # cl <- parallel::makeCluster(getOption("cl.cores", getCores()))
    step <- 50 # Avoid updating after analysing each event
    startProgress("Performing statistical analysis", 
                  divisions=5 + round(nrow(psi)/step))
    time <- Sys.time()
    
    count <- 0
    stats <- apply(psi, 1, function(...) {
        count <<- count + 1
        if (count %% step == 0)
            updateProgress("Performing statistical analysis", console=FALSE)
        return(statsAnalyses(...))
    }, factor(type), threshold=1, step=step, analyses=statsChoices)
    print(Sys.time() - time)
    
    # Check the column names of the different columns
    ns <- lapply(stats, names)
    uniq <- unique(ns)
    match <- fmatch(ns, uniq)
    
    # Convert matrix that share the same name to data frame (way faster)
    updateProgress("Preparing data")
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
        updateProgress("Calculating delta variance and median")
        deltaVar <- deltaVar[, 2] - deltaVar[, 1]
        deltaMed <- df[, grepl("Median", colnames(df))]
        deltaMed <- deltaMed[, 2] - deltaMed[, 1]
        df <- cbind(df, deltaVar, deltaMed)
    }
    
    if (any("density" == statsChoices)) {
        updateProgress("Calculating the density of inclusion levels")
        df[, "Density"] <- calculateDensitySparklines(df[, "Density"])
    }
    
    setDifferentialAnalyses(df)
    
    # parallel::stopCluster(cl)
    print(Sys.time() - time)
    closeProgress()
}

#' Optimal survival difference given an inclusion level cut-off for a specific
#' alternative splicing event
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
optimSurvDiff <- function(session, input, output) {
    observe({
        if (is.null(getDifferentialAnalyses()) || is.null(getClinicalData()))
            return(NULL)
        
        output$survivalOptions <- renderUI(optimSurvDiffUI(session$ns))
        
        # Update selectize input label depending on the chosen censoring type
        observe({
            if (is.null(input$censoring)) return(NULL)
            
            label <- "Follow up time"
            if (grepl("interval", input$censoring, fixed=TRUE))
                label <- "Starting time"
            updateSelectizeInput(session, "timeStart", label=label)
        })
        
        updateClinicalParams(session)
    })
    
    #' Calculate optimal survival cut-off for the inclusion levels of a given
    #' alternative splicing event
    observeEvent(input$survival, {
        survTerms <- isolate({
            # Get tumour sample IDs (normal and control samples are not
            # interesting for survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- getSampleTypes(names(match))
            tumour <- match[!grepl("Normal|Control", types)]

            # Group samples by the inclusion levels cut-off
            clinical <- getClinicalData()
            clinicalIDs <- nrow(clinical)
            groups <- rep(NA, clinicalIDs)

            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            dataEvent <- input$event
            psi       <- getInclusionLevels()
            stats     <- getDifferentialAnalyses()
            display   <- input$statsTable_rows_current
            selected  <- input$selected
        })

        if (selected == "shown") psi <- psi[display, ]
        startProgress("Performing survival analysis", nrow(psi))

        opt <- apply(psi, 1, function(vector) {
            v <- as.numeric(vector[toupper(names(tumour))])

            opt <- suppressWarnings(
                optim(0, testSurvivalCutoff, data=v, group=groups,
                      filter=tumour, clinical=clinical, censoring=censoring,
                      timeStart=timeStart, timeStop=timeStop, 
                      dataEvent=dataEvent, modals=FALSE,
                      # Method and parameters interval
                      method="Brent", lower=0, upper=1))

            updateProgress("Survival analysis", console=FALSE)
            return(c("Optimal survival PSI cut-off"=opt$par,
                     "Optimal survival difference"=opt$value))
        })
        # Remove NAs and add information to the statistical table
        opt <- opt[ , !is.na(opt[2, ])]
        df <- data.frame(t(opt))
        for (col in names(df)) stats[rownames(df), col] <- df[ , col]

        setDifferentialAnalyses(stats)
        closeProgress()
    })
}

#' Server logic of the exploratory differential analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom DT renderDataTable
diffAnalysisTableServer <- function(input, output, session) {
    ns <- session$ns
    
    # Information on the data groups from TCGA
    output$groupsInfo <- renderUI({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        
        if (is.null(psi)) return(tagList(
            helpText(icon("exclamation-circle"), 
                     "No alternative splicing junction quantification loaded.")))
        
        # Separate samples by their type
        ids <- names(psi)
        type <- getSampleTypes(ids)
        
        bullet <- "\u2022"
        groups <- NULL
        for (each in unique(type)) groups <- tagList(groups, br(), bullet, each)
        
        return(tagList(
            helpText("The data contains the following groups:", groups),
            hr()))
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            statsChoices <- input$statsChoices
        })
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        }
        
        performStatsAnalyses(session, input, psi, statsChoices)
    })
    
    # Show table and respective interface when statistical table is available
    observe({
        stats <- getDifferentialAnalyses()
        if (is.null(stats)) return(NULL)
        
        # Columns to show in statistical table
        output$showColumns <- renderUI({
            tagList(
                downloadButton(ns("download"), "Download whole table"),
                selectizeInput(ns("columns"), "Show columns", multiple=TRUE,
                               choices=colnames(stats), width="auto",
                               selected=colnames(stats),
                               options=list(plugins=list('remove_button', 
                                                         'drag_drop'))),
                hr())
        })
        
        # Render statistical table with the selected columns
        output$statsTable <- renderDataTableSparklines({
            if (!is.null(input$columns) && all(input$columns %in% names(stats)))
                stats[, input$columns]
        }, style="bootstrap", selection="none", filter='top',
        options=list(pageLength=10, rowCallback=JS("createDiffSplicingLinks")))
        
        # Prepare table to be downloaded
        output$download <- downloadHandler(
            filename=paste(getCategories(), "Differential splicing analyses"),
            content=function(file)
                write.table(getDifferentialAnalyses()[-c("Density")], 
                            file, quote=FALSE, row.names=FALSE, sep="\t")
        )
    })
    
    # Optimal survival difference given an inclusion level cut-off for a 
    # specific alternative splicing event
    optimSurvDiff(session, input, output)
}

attr(diffAnalysisTableUI, "loader") <- "analysis"
attr(diffAnalysisTableUI, "name") <- "Differential analysis (exploratory)"
attr(diffAnalysisTableServer, "loader") <- "analysis"