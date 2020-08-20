#' @include app.R
NULL

#' Missing information modal template
#'
#' @param session Shiny session
#' @param dataType Character: type of data missing
#' @param buttonId Character: identifier of button to take user to load missing
#' data
#'
#' @inherit psichomics return
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' if (shiny::isRunning()) {
#'     session <- session$ns
#'     buttonInput <- "takeMeThere"
#'     buttonId <- ns(buttonInput)
#'     dataType <- "Inclusion levels"
#'     missingDataModal(session, buttonId, dataType)
#'     observeEvent(input[[buttonInput]], missingDataGuide(dataType))
#' }
#' }
missingDataModal <- function(session, dataType, buttonId) {
    template <- function(buttonLabel) {
        errorModal(
            session, paste("Load", tolower(dataType)),
            "This analysis requires", tolower(dataType), "to proceed.",
            footer=actionButton(buttonId, buttonLabel, "data-dismiss"="modal",
                                class="btn-danger"),
            caller="Data analysis")
    }

    switch(dataType,
           "Inclusion levels"=template("Load or calculate"),
           template("Load"))
}

#' @rdname missingDataModal
missingDataGuide <- function(dataType) {
    js <- loadRequiredData(dataType)
    runjs(js)
}

#' Create input to select a gene
#'
#' @param id Character: identifier
#' @inheritParams shiny::selectizeInput
#' @param ... Arguments passed to the \code{options} list of
#'   \code{selectizeInput()}
#' @param placeholder Character: placeholder
#'
#' @return HTML elements
#' @keywords internal
selectizeGeneInput <- function(id, label="Gene", choices=NULL, multiple=FALSE,
                               ...,
                               placeholder="Type to search for a gene...") {
    onFocus  <- NULL
    onChange <- NULL
    if (!multiple) {
        onFocus  <- I('function() { this.clear(); }')
        onChange <- I('function(value) { this.blur(); }')
    }
    selectizeInput(
        id, label, width="100%", multiple=multiple, choices=choices,
        options=list(placeholder=placeholder,
                     plugins=list('remove_button', 'drag_drop'),
                     onFocus=onFocus, onChange=onChange, ...))
}

#' @rdname appUI
#'
#' @param id Character: identifier
#' @param tab Function to process HTML elements
#'
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS div icon fluidRow column tags
#'
#' @keywords internal
analysesUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(
        ns, "analysis",
        priority=paste0(c("dimReduction", "diffSplicing", "diffExpression",
                          "correlation", "survival", "info"), "UI"))

    # Load available analyses
    analyses <- tagList()
    for (k in seq(uiList)) {
        ui <- uiList[[k]]
        if ( is(ui, "shiny.tag.list") ) {
            tabName  <- attr(ui, "name")
            analyses <- c(analyses, list(tabPanel(tabName, ui)))
        } else {
            headerName <- attr(ui, "name")
            analyses   <- c(analyses, list(headerName))
            for (each in seq(ui)) {
                name     <- attr(ui[[each]], "name")
                analyses <- c(analyses, list(tabPanel(name, ui[each])))
            }
        }
        analyses <- c(analyses, list("----")) # Add a divider after each section
    }
    analyses <- analyses[-length(analyses)]
    do.call(tab, c(list(icon="flask", title="Analyses", menu=TRUE), analyses))
}

#' @rdname appServer
#'
#' @importFrom shiny observe observeEvent
#' @importFrom shinyjs hide show
#'
#' @keywords internal
analysesServer <- function(input, output, session) {
    # Run server logic from the scripts
    getServerFunctions("analysis",
                       priority=paste0(c("dimReduction", "diffSplicing",
                                         "diffExpression", "correlation",
                                         "survival", "info"),
                                       "Server"))
}

# Survival analyses helper functions --------------------------------------

#' Helper text to explain what happens when a subject matches multiple samples
#' when performing survival analysis
#'
#' @return Character
#' @keywords internal
subjectMultiMatchWarning <- function() {
    paste("While stratifying subjects for survival analysis, subjects",
          "with multipe samples are assigned the average value of their",
          "corresponding samples. However, for subjects with both disease",
          "and normal samples, it may be inappropriate to include the",
          "values of their normal samples for survival analysis.")
}

#' Retrieve clinical data based on attributes required for survival analysis
#'
#' @param ... Character: names of columns to retrieve
#' @param formulaStr Character: right-side of the formula for survival analysis
#'
#' @return Filtered clinical data
#' @keywords internal
getClinicalDataForSurvival <- function(..., formulaStr=NULL) {
    cols <- unlist(list(...))
    if (!is.null(formulaStr) && formulaStr != "") {
        form <- formula(paste("1 ~", formulaStr))
        cols <- c(cols, all.vars(form))
    }
    clinical <- getClinicalData(cols)
    return(clinical)
}

#' Assign average sample values to their corresponding subjects
#'
#' @param data One-row data frame/matrix or vector: values per sample for a
#' single gene
#' @param match Matrix: match between samples and subjects
#' @param clinical Data frame or matrix: clinical dataset (only required if the
#' \code{subjects} argument is not handed)
#' @param patients Character: subject identifiers (only required if the
#' \code{clinical} argument is not handed)
#' @param samples Character: samples to use when assigning values per subject
#' (if \code{NULL}, all samples will be used)
#'
#' @aliases getValuePerSubject getValuePerPatient assignValuePerPatient
#' @family functions to analyse survival
#' @return Values per subject
#' @export
#'
#' @examples
#' # Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#'
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#'
#' # Match between subjects and samples
#' match <- rep(paste("Subject", 1:3), 2)
#' names(match) <- colnames(psi)
#'
#' # Assign PSI values to each subject based on the PSI of their samples
#' assignValuePerSubject(psi[3, ], match)
assignValuePerSubject <- function(data, match, clinical=NULL, patients=NULL,
                               samples=NULL) {
    hasOneRow     <- !is.null(nrow(data)) && nrow(data) == 1
    isNamedVector <- is.vector(data) && !is.null(names(data))
    if (!hasOneRow && !isNamedVector)
        stop("Data needs to either have only one row or be a vector with ",
             "sample identifiers as names.")

    # TO DO: filter by subjects (allow to input no subjects, i.e. no filtering)
    if (is.null(patients)) patients <- rownames(clinical)

    if (!is.numeric(data)) {
        ns   <- names(data)
        data <- as.numeric(data)
        names(data) <- ns
    }

    # Filter by samples to use
    if (!is.null(samples)) match <- match[names(match) %in% samples]
    match <- match[!is.na(match)]

    # For each subject, assign the average value of its respective samples
    res <- sapply(split(data[names(match)], match), mean, na.rm=TRUE)
    return(res)
}

#' @export
getValuePerPatient <- assignValuePerSubject

#' @export
assignValuePerPatient <- assignValuePerSubject

#' @export
getValuePerSubject <- assignValuePerSubject

#' Process survival data to calculate survival curves
#'
#' @inheritParams getAttributesTime
#' @param group Character: group relative to each subject
#' @param clinical Data frame: clinical data
#' @param survTime \code{survTime} object: Times to follow up, time start, time
#' stop and event (optional)
#'
#' @inherit processSurvTerms details
#'
#' @return Data frame with terms needed to calculate survival curves
#' @keywords internal
processSurvData <- function(event, timeStart, timeStop, followup, group,
                            clinical, survTime=NULL) {
    if ( is.null(survTime) ) {
        survTime <- getAttributesTime(clinical, event, timeStart, timeStop,
                                      followup)
    }

    # Create new time using the starting time replacing the NAs with
    # days to last follow up
    nas <- is.na(survTime$start)
    survTime$time <- survTime$start
    survTime$time[nas] <- survTime$followup[nas]

    # Indicate event of interest and groups
    survTime$event <- ifelse(!is.na(survTime$event), 1, 0)
    if (!is.null(names(group)))
        survTime[names(group), "groups"] <- group
    else
        survTime$groups <- group

    if (!is.null(timeStop)) {
        # Create new time using the ending time replacing the NAs
        # with days to last follow up
        nas <- is.na(survTime$stop)
        survTime$time2 <- survTime$stop
        survTime$time2[nas] <- survTime$followup[nas]
    }
    return(survTime)
}

#' Get time values for given columns in a clinical dataset
#'
#' @param clinical Data frame: clinical data
#' @param event Character: name of column containing time of the event of
#' interest
#' @param timeStart Character: name of column containing starting time of the
#' interval or follow up time
#' @param timeStop Character: name of column containing ending time of the
#' interval (only relevant for interval censoring)
#' @param followup Character: name of column containing follow up time
#'
#' @family functions to analyse survival
#' @return Data frame containing the time for the given columns
#' @export
#'
#' @examples
#' df <- data.frame(followup=c(200, 300, 400), death=c(NA, 300, NA))
#' rownames(df) <- paste("subject", 1:3)
#' getAttributesTime(df, event="death", timeStart="death", followup="followup")
getAttributesTime <- function(clinical, event, timeStart, timeStop=NULL,
                              followup="days_to_last_followup") {
    cols <- c(followup=followup, start=timeStart, stop=timeStop, event=event)

    # Retrive time for given attributes
    timePerSubject <- function(col, clinical) {
        cols <- grep(col, colnames(clinical), value=TRUE)
        row  <- apply(clinical[cols], 1, function(i)
            if(!all(is.na(i))) max(as.numeric(i), na.rm = TRUE) else NA)
        return(row)
    }
    survTime <- lapply(cols, timePerSubject, clinical)

    survTime <- as.data.frame(survTime)
    class(survTime) <- c("data.frame", "survTime")
    return(survTime)
}

#' Update available clinical attributes when the clinical data changes
#'
#' @param session Shiny session
#' @param attrs Character: subject attributes
#'
#' @importFrom shiny observe updateSelectizeInput
#' @inherit psichomics return
#' @keywords internal
updateClinicalParams <- function(session, attrs) {
    if (!is.null(attrs)) {
        # Allow the user to select any "days_to" attribute available
        daysTo <- grep("days_to_", attrs, value=TRUE, fixed=TRUE)
        subDaysTo <- gsub(".*(days_to_.*)", "\\1", daysTo)
        choices <- unique(subDaysTo)
        names(choices) <- gsub("_", " ", choices, fixed=TRUE)
        names(choices) <- capitalize(names(choices))

        # Update choices for starting or follow up time
        updateSelectizeInput(
            session, "timeStart", choices=list(
                "Select a clinical attribute"="",
                "Suggested attributes"=choices,
                "All clinical attributes"=attrs),
            selected="days_to_death")

        # Update choices for ending time
        updateSelectizeInput(
            session, "timeStop", choices=list(
                "Select a clinical attribute"="",
                "Suggested attributes"=choices,
                "All clinical attributes"=attrs))

        # Update choices for events of interest
        names(choices) <- gsub("Days to ", "", names(choices), fixed=TRUE)
        names(choices) <- capitalize(names(choices))
        updateSelectizeInput(
            session, "event", choices=list(
                "Select a clinical attribute"="",
                "Suggested events"=choices,
                "All clinical attributes"=attrs),
            selected="days_to_death")
    } else {
        choices <- c("No clinical data loaded"="")
        updateSelectizeInput(session, "timeStart", choices=choices)
        updateSelectizeInput(session, "timeStop",  choices=choices)
        updateSelectizeInput(session, "event",     choices=choices)
    }
}

#' Process survival curves terms to calculate survival curves
#'
#' @inheritParams processSurvData
#' @param censoring Character: censor using \code{left}, \code{right},
#' \code{interval} or \code{interval2}
#' @param scale Character: rescale the survival time to \code{days},
#' \code{weeks}, \code{months} or \code{years}
#' @param formulaStr Character: formula to use
#' @param coxph Boolean: fit a Cox proportional hazards regression model?
#' @param survTime \code{survTime} object: times to follow up, time start, time
#' stop and event (optional)
#'
#' @importFrom stats formula
#' @importFrom survival coxph Surv
#'
#' @details The \code{event} time is only used to determine whether the event
#' has occurred (\code{1}) or not (\code{0}) in case of missing values.
#'
#' If \code{survTime = NULL}, survival times are obtained from the clinical
#' dataset according to the names given in \code{timeStart}, \code{timeStop},
#' \code{event} and \code{followup}. This may become quite slow when used in a
#' loop. If the aforementioned variables are constant, consider running
#' \code{\link{getAttributesTime}()} outside the loop and using its output via
#' the \code{survTime} argument of this function (see Examples).
#'
#' @family functions to analyse survival
#' @return A list with a \code{formula} object and a data frame with terms
#' needed to calculate survival curves
#' @export
#'
#' @examples
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup",
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#' formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
#' survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
#'                                formulaStr=formulaStr)
#'
#' # If running multiple times, consider calculating survTime only once
#' survTime <- getAttributesTime(clinical, event, timeStart)
#' for (i in seq(5)) {
#'   survTerms <- processSurvTerms(clinical, censoring="right", event,
#'                                 timeStart, formulaStr=formulaStr,
#'                                 survTime=survTime)
#' }
processSurvTerms <- function(clinical, censoring, event, timeStart,
                             timeStop=NULL, group=NULL, formulaStr=NULL,
                             coxph=FALSE, scale="days",
                             followup="days_to_last_followup", survTime=NULL) {
    # Ignore timeStop if interval-censoring is not selected
    if (!grepl("interval", censoring, fixed=TRUE) || timeStop == "")
        timeStop <- NULL

    # Check if using or not interval-censored data
    formulaSurv <- ifelse(is.null(timeStop),
                          "Surv(time/%s, event, type=censoring) ~",
                          "Surv(time/%s, time2, event, type=censoring) ~")
    scaleStr <- scale
    scale <- switch(scaleStr, days=1, weeks=7, months=30.42, years=365.25)
    formulaSurv <- sprintf(formulaSurv, scale)

    survData <- processSurvData(event, timeStart, timeStop, followup, group,
                                clinical, survTime)

    # Estimate survival curves by groups or using formula
    if (formulaStr == "" || is.null(formulaStr)) {
        formulaTerms <- "groups"
    } else {
        formulaTerms <- formulaStr
        factor <- sapply(clinical, is.factor)
        clinical[factor] <- lapply(clinical[factor], as.character)
        survData <- cbind(survData, clinical)
    }

    form <- formula(paste(formulaSurv, formulaTerms))

    if (coxph) {
        res <- coxph(form, data=survData)

        if (!is.null(res$xlevels$groups)) {
            # Correct group names
            name <- res$xlevels$groups[-1]
            names(res$means) <- name
            names(res$coefficients)  <- name
            # names(res$wald.test)  <- name
        }
    } else {
        res <- list(form=form, survTime=survData)
    }

    res$scale <- scaleStr
    class(res) <- c("survTerms", class(res))
    return(res)
}

#' @inherit survival::survfit title details
#' @inheritParams survdiffTerms
#' @inheritDotParams survival::survdiff -formula -data
#'
#' @importFrom survival survfit
#'
#' @family functions to analyse survival
#' @inherit survdiffTerms return
#' @export
#'
#' @examples
#' library("survival")
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup",
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#' formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
#' survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
#'                                formulaStr=formulaStr)
#' survfit(survTerms)
survfit.survTerms <- function(survTerms, ...) {
    res <- survfit(survTerms$form, data=survTerms$survTime, ...)
    res$scale <- survTerms$scale

    # Correct group names
    groups <- deparse(survTerms$form[[3]])
    if (!is.null(res$strata) && groups == "groups") {
        name <- paste0("^", groups, "=")
        names(res$strata) <- gsub(name, "", names(res$strata))
    }
    return(res)
}

#' @inherit survival::survdiff
#' @param survTerms \code{survTerms} object: survival terms obtained after
#'   running \code{processSurvTerms} (see examples)
#' @inheritDotParams survival::survdiff -formula -data
#'
#' @importFrom survival survdiff
#'
#' @family functions to analyse survival
#' @return \code{survfit} object. See \code{survfit.object} for details. Methods
#' defined for \code{survfit} objects are \code{print}, \code{plot},
#' \code{lines}, and \code{points}.
#' @export
#'
#' @examples
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup",
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#' formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
#' survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
#'                                formulaStr=formulaStr)
#' survdiffTerms(survTerms)
survdiffTerms <- function(survTerms, ...) {
    survdiff(survTerms$form, data=survTerms$survTime, ...)
}

#' Plot survival curves
#'
#' @param surv Survival object
#' @param interval Boolean: show interval ranges?
#' @param mark Boolean: mark times?
#' @param title Character: plot title
#' @param pvalue Numeric: p-value of the survival curves
#' @param scale Character: time scale (default is \code{days})
#' @param auto Boolean: return the plot automatically prepared (\code{TRUE}) or
#' only the bare minimum (\code{FALSE})?
#'
#' @importFrom shiny tags br
#'
#' @family functions to analyse survival
#' @return Plot of survival curves
#' @export
#'
#' @examples
#' require("survival")
#' fit <- survfit(Surv(time, status) ~ x, data = aml)
#' plotSurvivalCurves(fit)
plotSurvivalCurves <- function(surv, mark=TRUE, interval=FALSE, pvalue=NULL,
                               title="Survival analysis", scale=NULL,
                               auto=TRUE) {
    hc <- hchart(surv, ranges=interval, markTimes=mark)
    if (auto) {
        if (is.null(scale)) {
            if (is.null(surv$scale))
                scale <- "days"
            else
                scale <- surv$scale
        }

        hc <- hc %>%
            hc_chart(zoomType="xy") %>%
            hc_title(text=unname(title)) %>%
            hc_yAxis(title=list(text="Survival proportion"),
                     crosshair=TRUE) %>%
            hc_xAxis(title=list(text=paste("Time in", scale)),
                     crosshair=TRUE) %>%
            hc_tooltip(
                headerFormat = paste(
                    tags$small("{point.x:.3f}", scale), br(),
                    span(style="color:{point.color}", "\u25CF "),
                    tags$b("{series.name}"), br()),
                pointFormat = paste(
                    "Survival proportion: {point.y:.3f}", br(),
                    "Records: {series.options.records}", br(),
                    "Events: {series.options.events}", br(),
                    "Median: {series.options.median}")) %>%
            hc_plotOptions(series=list(stickyTracking=FALSE))
    }

    if (!is.null(pvalue))
        hc <- hc_subtitle(hc, text=paste("log-rank p-value:", pvalue))
    return(hc)
}

#' Check if survival analyses successfully completed or returned errors
#'
#' @param session Shiny session
#' @inheritDotParams processSurvTerms
#'
#' @importFrom shiny tags
#'
#' @return List with survival analysis results
#' @keywords internal
processSurvival <- function(session, ...) {
    # Calculate survival curves
    survTerms <- tryCatch(processSurvTerms(...), error=return)
    if ("simpleError" %in% class(survTerms)) {
        if (survTerms[[1]] == paste("contrasts can be applied only to",
                                    "factors with 2 or more levels")) {
            errorModal(session, "Formula error",
                       "Cox models can only be applied to 2 or more groups.",
                       caller="Data analysis")
        } else {
            errorModal(session, "Formula error",
                       "Maybe you misplaced a ", tags$kbd("+"), ", ",
                       tags$kbd(":"), " or ", tags$kbd("*"), "?", br(),
                       br(), "The following error was raised:", br(),
                       tags$code(survTerms$message),
                       caller="Data analysis")
        }
        return(NULL)
    }
    return(survTerms)
}

#' Test the survival difference between groups of subjects
#'
#' @inheritParams survdiffTerms
#' @inheritDotParams survival::survdiff -formula -data
#'
#' @note Instead of raising errors, returns \code{NA}
#'
#' @family functions to analyse survival
#' @return p-value of the survival difference or \code{NA}
#' @export
#'
#' @examples
#' require("survival")
#' data <- aml
#' timeStart  <- "event"
#' event      <- "event"
#' followup   <- "time"
#' data$event <- NA
#' data$event[aml$status == 1] <- aml$time[aml$status == 1]
#' censoring  <- "right"
#' formulaStr <- "x"
#' survTerms <- processSurvTerms(data, censoring=censoring, event=event,
#'                               timeStart=timeStart, followup=followup,
#'                               formulaStr=formulaStr)
#' testSurvival(survTerms)
testSurvival <- function (survTerms, ...) {
    # If there's an error with survdiff, return NA
    pvalue <- tryCatch({
        # Test the difference between survival curves
        diff <- survdiffTerms(survTerms, ...)

        # Calculate p-value with given significant digits
        pvalue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
        return(as.numeric(signifDigits(pvalue)))
    }, error = function(e) NA)
    return(pvalue)
}

#' Label groups based on a given cutoff
#'
#' @param data Numeric: test data
#' @param cutoff Numeric: test cutoff
#' @param label Character: label to prefix group names
#' @param gte Boolean: test using greater than or equal than cutoff
#' (\code{TRUE}) or less than or equal than cutoff (\code{FALSE})?
#'
#' @family functions to analyse survival
#' @return Labelled groups
#' @export
#'
#' @examples
#' labelBasedOnCutoff(data=c(1, 0, 0, 1, 0, 1), cutoff=0.5)
#'
#' labelBasedOnCutoff(data=c(1, 0, 0, 1, 0, 1), cutoff=0.5, "Ratio")
#'
#' # Use "greater than" instead of "greater than or equal to"
#' labelBasedOnCutoff(data=c(1, 0, 0, 0.5, 0, 1), cutoff=0.5, gte=FALSE)
labelBasedOnCutoff <- function (data, cutoff, label=NULL, gte=TRUE) {
    len <- length(data)
    group <- rep(NA, len)

    if (gte) {
        comp <- `>=`
        str1 <- "&gt;="
        str2 <- "&lt;"
    } else {
        comp <- `>`
        str1 <- "&gt;"
        str2 <- "&lt;="
    }
    group <- comp(data, cutoff)

    # Assign a value based on the inclusion levels cutoff
    if (is.null(label)) {
        group[group == "TRUE"]  <- paste(str1, cutoff)
        group[group == "FALSE"] <- paste(str2, cutoff)
    } else {
        group[group == "TRUE"]  <- paste(label, str1, cutoff)
        group[group == "FALSE"] <- paste(label, str2, cutoff)
    }

    length(group) <- len
    return(group)
}

#' Test the survival difference between two survival groups given a cutoff
#'
#' @inheritParams processSurvTerms
#' @param cutoff Numeric: Cutoff of interest
#' @param data Numeric: elements of interest to test against the cutoff
#' @param filter Boolean or numeric: elements to use (all are used by default)
#' @inheritDotParams processSurvTerms -group -clinical
#' @param session Shiny session
#' @param survivalInfo Boolean: return extra survival information
#'
#' @importFrom survival survdiff
#' @return p-value of the survival difference
#' @keywords internal
testSurvivalCutoff <- function(cutoff, data, filter=TRUE, clinical, ...,
                               session=NULL, survivalInfo=FALSE) {
    group <- labelBasedOnCutoff(data, cutoff, label="Inclusion levels")

    # Calculate survival curves
    if (!is.null(session)) {
        survTerms <- processSurvival(session, group=group, clinical=clinical,
                                     ...)
        if (is.null(survTerms)) return(NULL)
    } else {
        survTerms <- tryCatch(processSurvTerms(group=group, clinical=clinical,
                                               ...),
                              error=return)
        if ("simpleError" %in% class(survTerms)) return(NA)
    }

    pvalue <- testSurvival(survTerms)
    if (is.na(pvalue)) pvalue <- 1
    if (survivalInfo) attr(pvalue, "info") <- survfit(survTerms)
    return(pvalue)
}

#' Calculate optimal data cutoff that best separates survival curves
#'
#' Uses \code{stats::optim} with the Brent method to test multiple cutoffs and
#' to find the minimum log-rank p-value.
#'
#' @inheritParams processSurvTerms
#' @inheritParams testSurvivalCutoff
#' @param data Numeric: data values
#' @param session Shiny session (only used for the visual interface)
#' @param lower,upper Bounds in which to search (if \code{NULL}, bounds are set
#' to \code{lower = 0} and \code{upper = 1} if all data values are within that
#' interval; otherwise, \code{lower = min(data, na.rm = TRUE)} and
#' \code{upper = max(data, na.rm = TRUE)})
#'
#' @family functions to analyse survival
#' @return List containing the optimal cutoff (\code{par}) and the corresponding
#' p-value (\code{value})
#' @export
#'
#' @examples
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup",
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#'
#' psi <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
#' opt <- optimalSurvivalCutoff(clinical, psi, "right", event, timeStart)
optimalSurvivalCutoff <- function(clinical, data, censoring, event, timeStart,
                                  timeStop=NULL,
                                  followup="days_to_last_followup",
                                  session=NULL, filter=TRUE, survTime=NULL,
                                  lower=NULL, upper=NULL) {
    if (is.null(lower) && is.null(upper)) {
        # Search between min and max of data
        lower <- min(data, na.rm=TRUE)
        upper <- max(data, na.rm=TRUE)

        if (lower >= 0 && upper <= 1) {
            # Search between 0 and 1 (if data values are within that interval)
            lower <- 0
            upper <- 1
        } else if (lower >= upper) {
            upper <- lower + 1
        }
    }

    if ( is.null(survTime) )
        survTime <- getAttributesTime(clinical, event, timeStart, timeStop,
                                      followup)

    # Supress warnings from failed calculations while optimising
    opt <- suppressWarnings(
        optim(0, testSurvivalCutoff, data=data, filter=filter,
              clinical=clinical, censoring=censoring, timeStart=timeStart,
              timeStop=timeStop, event=event, followup=followup,
              survTime=survTime, session=session,
              # Method and parameters interval
              method="Brent", lower=lower, upper=upper))
    return(opt)
}

# Differential analyses helper functions -----------------------------------

#' Prepare event plot options
#'
#' @param id Character: identifier
#' @param ns Namespace identifier
#' @param labelsPanel Tab panel containing options to label points
#'
#' @return HTML elements
#' @keywords internal
prepareEventPlotOptions <- function(id, ns, labelsPanel=NULL) {
    createAxisPanel <- function(ns, axis) {
        upper  <- toupper(axis)
        idAxis <- function(label) ns(paste0(axis, label))
        highlightLabel <- sprintf("Highlight points based on %s values", upper)

        tab <- tabPanel(
            paste(upper, "axis"),
            selectizeInput(idAxis("Axis"), choices=NULL, width="100%",
                           sprintf("Select %s axis", upper)),
            selectizeInput(idAxis("Transform"),
                           sprintf("Data transformation of %s values", upper),
                           transformOptions(axis), width="100%"),
            bsCollapse(bsCollapsePanel(
                list(icon("thumb-tack"), highlightLabel),
                value=paste0(axis, "AxisHighlightPanel"),
                checkboxInput(idAxis("Highlight"), width="100%",
                              label=highlightLabel),
                uiOutput(idAxis("HighlightValues")))))
        return(tab)
    }

    xAxisPanel <- createAxisPanel(ns, "x")
    yAxisPanel <- createAxisPanel(ns, "y")

    plotStyle <- navbarMenu(
        "Plot style",
        tabPanel("Base points",
                 plotPointsStyle(
                     ns, "base", "Base points",
                     help="These are points not highlighted or selected.",
                     size=2, colour="gray", alpha=0.3)),
        tabPanel("Highlighted points",
                 plotPointsStyle(
                     ns, "highlighted", "Highlighted points",
                     help="Highlight points in the X and Y axes options.",
                     size=3, colour="orange", alpha=0.5)),
        tabPanel("Selected in the table",
                 plotPointsStyle(
                     ns, "selected", "Selected in the table",
                     help=paste("Click in a row of the table to emphasise the",
                                "respective point in the plot."),
                     size=8, colour="blue", alpha=0.5)),
        tabPanel("Labels",
                 plotPointsStyle(
                     ns, "labelled", "Labels",
                     help="Modify the style of labels",
                     size=4, colour="red", alpha=1)))

    div(id=id, tabsetPanel(xAxisPanel, yAxisPanel, labelsPanel, plotStyle))
}

#' Perform and display statistical analysis
#'
#' Includes interface containing the results
#'
#' @inheritParams plotDistribution
#' @param stat Data frame or matrix: values of the analyses to be performed (if
#' \code{NULL}, the analyses will be performed)
#'
#' @details
#' \itemize{
#'   \item{\code{ttest}: unpaired t-test}
#'   \item{\code{wilcox}: Wilcoxon test}
#'   \item{\code{levene}: Levene's test}
#'   \item{\code{fligner}: Fligner-Killeen test}
#'   \item{\code{kruskal}: Kruskal test}
#'   \item{\code{fisher}: Fisher's exact test}
#'   \item{\code{spearman}: Spearman's test}
#' }
#'
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats wilcox.test
#' @importFrom R.utils capitalize
#'
#' @return HTML elements
#' @keywords internal
wilcox <- function(data, groups, stat=NULL) {
    warn <- NULL
    group <- unique(groups)
    len <- length(group)

    p.value <- NULL
    if (!is.null(stat)) {
        method      <- stat$`Wilcoxon method`
        statistic   <- stat$`Wilcoxon statistic`
        p.value     <- stat$`Wilcoxon p-value`
        null.value  <- stat$`Wilcoxon null value`
        alternative <- stat$`Wilcoxon alternative`
    }

    if (len != 2) {
        return(tagList(h4("Wilcoxon test"),
                       "Can only perform this test on 2 groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Wilcoxon p-value \\(.* adjusted\\)", colnames(stat),
                         value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        dataA <- data[groups == group[1]]
        dataB <- data[groups == group[2]]
        stat <- tryCatch(list(stat=wilcox.test(dataA, dataB)),
                         warning=function(w)
                             return(list(stat=wilcox.test(dataA, dataB),
                                         warning=w)))

        if ("warning" %in% names(stat))
            warn <- tags$div(class="alert alert-warning", role="alert",
                             capitalize(stat$warning$message))

        method      <- stat$stat$method
        statistic   <- stat$stat$statistic
        p.value     <- stat$stat$p.value
        adjusted    <- NULL
        null.value  <- stat$stat$null.value
        alternative <- stat$stat$alternative
    }

    tagList(
        h4(method), warn,
        tags$b("Test value: "), roundDigits(statistic), br(),
        tags$b("Location parameter: "), null.value, br(),
        tags$b("Alternative hypothesis: "), alternative,
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' @rdname wilcox
#'
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats t.test
#' @importFrom R.utils capitalize
ttest <- function(data, groups, stat=NULL) {
    warn <- NULL
    group <- unique(groups)
    len <- length(group)

    p.value <- NULL
    if (!is.null(stat)) {
        method      <- stat$`T-test method`
        statistic   <- stat$`T-test statistic`
        p.value     <- stat$`T-test p-value`
        null.value  <- stat$`T-test null value`
        alternative <- stat$`T-test alternative`
        parameter   <- stat$`T-test parameter`
        int1        <- stat$`T-test conf int1`
        int2        <- stat$`T-test conf int2`
    }

    if (len != 2) {
        return(tagList(h4("Unpaired t-test"),
                       "Can only perform this test on 2 groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("T-test p-value \\(.* adjusted\\)", colnames(stat),
                         value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        dataA <- data[groups == group[1]]
        dataB <- data[groups == group[2]]
        stat <- tryCatch(list(stat=t.test(dataA, dataB)),
                         warning=function(w)
                             return(list(stat=t.test(dataA, dataB),
                                         warning=w)),
                         error=return)
        if (is(stat, "error")) {
            message <- stat$message
            check <- "not enough '%s' observations"
            checkX <- sprintf(check, "x")
            checkY <- sprintf(check, "y")

            fewObservations <- function(name)
                tagList("Not enough observations in group", tags$b(name),
                        "to perform this statistical test.")

            if (message == checkX)
                message <- fewObservations(group[1])
            else if (message == checkY)
                message <- fewObservations(group[2])
            else
                message <- capitalize(message)

            error <- tagList(h4("t-test"), tags$div(class="alert alert-danger",
                                                    role="alert", message))
            return(error)
        }

        if ("warning" %in% names(stat))
            warn <- tags$div(class="alert alert-warning", role="alert",
                             capitalize(stat$warning$message))

        method      <- stat$stat$method
        statistic   <- stat$stat$statistic
        p.value     <- stat$stat$p.value
        adjusted    <- NULL
        null.value  <- stat$stat$null.value
        alternative <- stat$stat$alternative
        parameter   <- stat$stat$parameter
        int1        <- stat$stat$conf.int[[1]]
        int2        <- stat$stat$conf.int[[2]]
    }

    tagList(
        h4(method), warn,
        tags$b("Test value: "), roundDigits(statistic), br(),
        tags$b("Test parameter: "), parameter, br(),
        tags$b("Difference in means: "), null.value, br(),
        tags$b("Alternative hypothesis: "), alternative, br(),
        tags$b("95\u0025 confidence interval: "), roundDigits(int1),
        roundDigits(int2),
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' @rdname wilcox
#' @importFrom shiny tagList tags h4 br
levene <- function(data, groups, stat=NULL) {
    p.value <- NULL
    if (!is.null(stat)) {
        statistic <- stat$`Levene statistic`
        p.value   <- stat$`Levene p-value`
        non.bootstrap.p.value <- stat$`Levene non bootstrap p-value`
    }

    len <- length(unique(groups))
    if (len < 2) {
        return(tagList(h4("Levene's Test for Homogeneity of Variances"),
                       "Can only perform this test on 2 or more groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Levene .*p-value \\(.* adjusted\\)", colnames(stat),
                         value=TRUE)
        if (length(adjusted) == 2) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label1 <- paste0("p-value (", adjustMethod[[1]], "): ")
            label2 <- paste0("p-value without bootstrap (", adjustMethod[[2]],
                             "): ")
            adjustedNonBootstrap <- tagList(br(), tags$b(label2),
                                            signifDigits(adjusted[[1]]), br())
            adjusted <- tagList(tags$b(label1), signifDigits(adjusted[[2]]),
                                br())
        } else if (length(adjusted) == 1) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
            adjustedNonBootstrap <- NULL
        }
    } else {
        nas <- is.na(data)
        stat <- leveneTest(data[!nas], factor(groups[!nas]))
        statistic <- stat$statistic
        p.value   <- stat$p.value
        adjusted  <- NULL
        non.bootstrap.p.value <- stat$non.bootstrap.p.value
        adjustedNonBootstrap  <- NULL
    }

    if (!is.null(non.bootstrap.p.value)) {
        nonBootstrap <- tagList(tags$b("p-value without bootstrap: "),
                                signifDigits(non.bootstrap.p.value),
                                adjustedNonBootstrap)
    } else {
        nonBootstrap <- NULL
    }

    tagList(
        h4("Levene's Test for Homogeneity of Variance"),
        tags$b("Test value: "), roundDigits(statistic), br(),
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted,
            nonBootstrap))
}

#' @rdname wilcox
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats fligner.test
fligner <- function(data, groups, stat=NULL) {
    len <- length(unique(groups))

    p.value <- NULL
    if (!is.null(stat)) {
        statistic <- stat$`Fligner-Killeen statistic`
        p.value   <- stat$`Fligner-Killeen p-value`
        parameter <- stat$`Fligner-Killeen parameter`
    }

    if (len < 2) {
        return(tagList(h4("Fligner-Killeen Test for Homogeneity of Variances"),
                       "Can only perform this test on 2 or more groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Fligner-Killeen .*p-value \\(.* adjusted\\)",
                         colnames(stat), value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        nas <- is.na(data)
        stat <- fligner.test(data[!nas], factor(groups[!nas]))
        if (is.infinite(stat$statistic)) {
            statistic <- NaN
            p.value   <- NA
            parameter <- NaN
        } else {
            statistic <- stat$statistic
            p.value   <- stat$p.value
            parameter <- stat$parameter
        }
        adjusted <- NULL
    }

    tagList(
        h4("Fligner-Killeen's Test for Homogeneity of Variance"),
        tags$b("Test value: "), roundDigits(statistic), br(),
        tags$b("Test parameter: "), parameter, br(),
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' @rdname wilcox
#'
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats kruskal.test
kruskal <- function(data, groups, stat=NULL) {
    len <- length(unique(groups))

    p.value <- NULL
    if (!is.null(stat)) {
        method    <- stat$`Kruskal method`
        statistic <- stat$`Kruskal statistic`
        p.value   <- stat$`Kruskal p-value`
        parameter <- stat$`Kruskal parameter`
    }

    if (len < 2) {
        return(tagList(h4("Kruskal test"),
                       "Can only perform this test on 2 or more groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Kruskal p-value \\(.* adjusted\\)", colnames(stat),
                         value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        stat      <- kruskal.test(data, factor(groups))
        method    <- stat$method
        statistic <- stat$statistic
        p.value   <- stat$p.value
        parameter <- stat$parameter
        adjusted  <- NULL
    }

    tagList(h4(method),
            tags$b("Test value \u03C7\u00B2: "), roundDigits(statistic), br(),
            tags$b("Degrees of freedom: "), parameter,
            div(style="text-align:right",
                tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' @rdname wilcox
#'
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats fisher.test
#' @importFrom R.utils withTimeout
fisher <- function(data, groups) {
    stat <- try(withTimeout(
        fisher.test(data, factor(groups)),
        timeout = 1,
        onTimeout = "error"))

    if (!is(stat, "try-error")) {
        tagList(
            h4(stat$method),
            tags$b("p-value: "), stat$p.value, br(),
            tags$b("Alternative hypothesis: "), stat$alternative
        )
    } else {
        tagList(
            h4("Fisher's Exact Test for Count Data"),
            "This test took too much to complete!"
        )
    }
}

#' @rdname wilcox
#'
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats var cor
spearman <- function(data, groups) {
    group <- unique(groups)
    len <- length(group)

    if (len != 2) {
        tagList(
            h4("Spearman's correlation"),
            "Can only perform this test on 2 groups.")
    } else {
        var <- var(data[groups == group[1]], data[groups == group[2]])
        cor <- cor(data[groups == group[1]], data[groups == group[2]])

        tagList(
            h4("Spearman's correlation"),
            tags$b("Variance: "), var, br(),
            tags$b("Correlation: "), cor)
    }
}

#' Options for event plotting
#'
#' @param session Shiny session
#' @param df Data frame
#' @param xAxis Character: currently selected variable for the X axis
#' @param yAxis Character: currently selected variable for the Y axis
#' @param labelSortBy Character: currently selected variable for the
#' \code{selectize} element to sort differentially analysis
#'
#' @return HTML elements
#' @keywords internal
eventPlotOptions <- function(session, df, xAxis, yAxis, labelSortBy) {
    # Only allow to select numeric columns
    cols <- colnames(df)
    type <- sapply(cols, function(i) class(df[[i]]))
    numericCols <- cols[type == "numeric"]

    if (!is.null(numericCols) && length(numericCols) > 0) {
        if (is.null(xAxis) || identical(xAxis, "") || !xAxis %in% numericCols) {
            # Default option for X axis
            deltaMedian <- "\u2206 Median"    # PSI
            logFC       <- "log2 Fold-Change" # Gene expression
            if (logFC %in% numericCols)
                xSelected <- logFC
            else if (deltaMedian %in% numericCols)
                xSelected <- deltaMedian
            else
                xSelected <- NULL

            updateSelectizeInput(session, "xAxis", choices=numericCols,
                                 selected=xSelected)
        } else {
            updateSelectizeInput(session, "xAxis", choices=numericCols,
                                 selected=xAxis)
        }

        if (is.null(yAxis) || identical(yAxis, "") || !yAxis %in% numericCols) {
            # Default option for Y axis
            pValue <- grepl("p-value", numericCols)
            bStat  <- "B-statistics" # Gene expression
            if (bStat %in% numericCols)
                ySelected <- bStat
            else if (any(pValue))
                ySelected <- numericCols[pValue][1]
            else if (length(numericCols) >= 2)
                ySelected <- numericCols[[2]]
            else
                ySelected <- NULL

            updateSelectizeInput(session, "yAxis", choices=numericCols,
                                 selected=ySelected)
        } else {
            updateSelectizeInput(session, "yAxis", choices=numericCols,
                                 selected=yAxis)
        }

        if (is.null(labelSortBy) || identical(labelSortBy, "") ||
            !labelSortBy %in% numericCols) {
            # Default option for sorting differentially analysis
            pValue <- grepl("p-value", numericCols)
            bStat  <- "B-statistics" # Gene expression
            if (bStat %in% numericCols)
                labelSelected <- bStat
            else if (any(pValue))
                labelSelected <- numericCols[pValue][1]
            else if (length(numericCols) >= 2)
                labelSelected <- numericCols[[2]]
            else
                labelSelected <- NULL

            updateSelectizeInput(session, "labelSortBy", choices=numericCols,
                                 selected=labelSelected)
        } else {
            updateSelectizeInput(session, "labelSortBy", choices=numericCols,
                                 selected=labelSortBy)
        }
    } else {
        updateSelectizeInput(session, "xAxis", choices=character(0))
        updateSelectizeInput(session, "yAxis", choices=character(0))
        updateSelectizeInput(session, "labelSortBy", choices=character(0))
    }
}

#' Basic statistics performed on data
#'
#' Variance and median of each group. If data has 2 groups, also calculates the
#' delta variance and delta median.
#'
#' @inheritParams plotDistribution
#' @importFrom shiny tagList br h4
#'
#' @return HTML elements
#' @keywords internal
basicStats <- function(data, groups) {
    data <- lapply(unique(groups), function(g) data[groups == g])

    len <- length(unique(groups))
    vari <- vapply(data, var, numeric(1), na.rm = TRUE)
    medi <- vapply(data, median, numeric(1), na.rm = TRUE)

    if (len == 2) {
        deltaMedian <- tagList(tags$b("|\u0394 Median|: "),
                               roundDigits(abs(medi[2] - medi[1])), br())
        deltaVar <- tagList(tags$b("|\u0394 Variance|: "),
                            roundDigits(abs(vari[2] - vari[1])), br())
    } else {
        deltaMedian <- NULL
        deltaVar <- NULL
    }

    avgMedian <- roundDigits( mean(medi) )
    avgVar <- roundDigits( mean(vari) )
    ui <- tagList(hr(), h4("Basic statistics"),
                  tags$b("Average median: "), avgMedian, br(), deltaMedian,
                  tags$b("Average variance: "), avgVar, br(), deltaVar)
    return(ui)
}

#' Filter groups with less data points than the threshold
#'
#' Groups containing a number of non-missing values less than the threshold are
#' discarded.
#'
#' @param vector Unnamed elements
#' @param group Character: group of the elements
#' @param threshold Integer: number of valid non-missing values by group
#'
#' @return Named vector with filtered elements from valid groups. The group of
#' the respective element is given in the name.
#' @export
#'
#' @examples
#' # Removes groups with less than two elements
#' filterGroups(1:4, c("A", "B", "B", "D"), threshold=2)
filterGroups <- function(vector, group, threshold=1) {
    names(vector) <- group
    vector <- lapply(unique(group), function(t) {
        vector <- vector[group == t]
        if ( sum(!is.na(vector)) >= threshold )
            return(vector)
    })
    return(unlist(vector))
}

#' Create plot for events
#'
#' @param df Data frame
#' @param x Character: name of the variable used for the X axis
#' @param y Character: name of the variable used for the Y axis
#' @param params List of parameters to pass to
#' \code{\link[ggplot2]{geom_point}()} related to most points
#' @param highlightX Integer: region of points in X axis to highlight
#' @param highlightY Integer: region of points in Y axis to highlight
#' @param highlightParams List of parameters to pass to
#' \code{\link[ggplot2]{geom_point}()} related to highlighted points
#' @param selected Integer: index of rows/points to be coloured
#' @param selectedParams List of parameters to pass to
#' \code{\link[ggplot2]{geom_point}()} related to selected points
#' @param labelled Integer: index of rows/points to be labelled
#' @param labelledParams List of parameters to pass to
#' \code{ggrepel::geom_label_repel} related to labelled points
#' @param xlim Numeric: limits of X axis
#' @param ylim Numeric: limits of Y axis
#'
#' @importFrom ggplot2 ggplot aes_string geom_point theme_light coord_cartesian
#' unit
#' @importFrom ggrepel geom_label_repel
#'
#' @return List containing HTML elements and highlighted points
#' @keywords internal
createEventPlotting <- function(df, x, y, params, highlightX, highlightY,
                                highlightParams, selected, selectedParams,
                                labelled, labelledParams, xlim, ylim) {
    aes <- aes_string(paste0("`", x, "`"), paste0("`", y, "`"))

    # Get points highlighted in X and Y that were not selected
    getNonSelectedPoints <- function(highlight, df, axis) {
        if (!is.null(highlight)) {
            onlyUpLimit <- identical(names(highlight), "upper")
            highlighted <- findInterval(df[[axis]], highlight, left.open=FALSE,
                                        rightmost.closed=TRUE) == !onlyUpLimit
            if (attr(highlight, "inverted")) highlighted <- !highlighted
            highlighted <- which(highlighted)
        } else {
            highlighted <- seq(nrow(df))
        }
        return(highlighted)
    }

    highlightedX <- getNonSelectedPoints(highlightX, df, x)
    highlightedY <- getNonSelectedPoints(highlightY, df, y)

    if ( is.null(highlightX) && is.null(highlightY) ) {
        highlighted <- NULL
    } else {
        highlighted <- intersect(highlightedX, highlightedY)
    }

    # Render remaining points
    plotted <- union(selected, highlighted)
    if (!is.null(plotted)) {
        remaining <- df[-plotted, ]
    } else {
        remaining <- df
    }
    plot <- ggplot() + do.call("geom_point", c(
        list(data=remaining, aes, na.rm=TRUE), params))

    # Render highlighted points
    plot <- plot + do.call("geom_point", c(
        list(data=df[setdiff(highlighted, selected), ], aes, na.rm=TRUE),
        highlightParams))

    # Render selected points
    plot <- plot + do.call("geom_point", c(
        list(data=df[selected, ], aes, na.rm=TRUE), selectedParams))

    # Label points
    aesMod <- aes
    aesMod$label <- parse(text="`Names`")[[1]]

    modNames <- rownames(df)
    if ( isTRUE(attr(labelled, "displayOnlyGene")) )
        modNames <- gsub(".*_(.*)$", "\\1", modNames)

    mod <- cbind(Names=modNames, df)
    plot <- plot + do.call("geom_label_repel", c(
        list(data=mod[labelled, ], aesMod, na.rm=TRUE), labelledParams))

    plot <- plot + coord_cartesian(xlim=xlim, ylim=ylim) + theme_light(16)
    return(list(plot=list(plot), highlighted=highlighted))
}

#' Show variable transformation(s)
#'
#' @param label Character: label to display
#' @param type Character: show the variable transformation for the chosen type;
#' if \code{NULL}, show all variable transformations
#'
#' @return Character labelling variable transformation(s)
#' @keywords internal
transformOptions <- function(label, type=NULL) {
    transform <- c("No transformation"="no",
                   "|%s|"="abs",
                   "-%s"="inv",
                   "log10(|%s|)"="log10abs",
                   "-log10(|%s|)"="-log10abs")
    names(transform) <- sprintf(names(transform), label)

    if (!is.null(type)) {
        show <- names(transform)
        show[[1]] <- label
        show <- show[match(type, transform)]
        return(show)
    } else {
        return(transform)
    }
}

#' Transform values as per a given type of transformation
#'
#' @param val Integer: values to transform
#' @param type Character: type of transformation
#' @param avoidZero Boolean: add the smallest non-zero number available
#' (\code{.Machine$double.xmin}) to avoid infinity values following
#' log-transformation (may not be plotted); useful for p-values of 0
#'
#' @return Integer containing transformed values
#' @keywords internal
transformValues <- function(val, type, avoidZero=TRUE) {
    # Remove NAs
    if (avoidZero) {
        zeroes <- val == 0 & !is.na(val)
        val[zeroes] <- val[zeroes] + .Machine$double.xmin
    }

    trans <- suppressWarnings(
        switch(type,
               "no"=val,
               "abs"=abs(val),
               "inv"=-val,
               "log10abs"=log10(abs(val)),
               "-log10abs"=-log10(abs(val))))
    return(trans)
}

#' Transform data in data frame
#'
#' @inheritParams appServer
#' @param df Data frame
#' @param x Character: column name
#' @param y Character: column name
#'
#' @return Data frame with transformed data in new columns and respective name
#' of created columns
#' @keywords internal
transformData <- function(input, df, x, y) {
    xTrans <- input$xTransform
    xLabel <- transformOptions(x, xTrans)
    if (!x %in% colnames(df)) return(NULL)
    df[[xLabel]] <- transformValues(df[[x]], xTrans)

    yTrans <- input$yTransform
    yLabel <- transformOptions(y, yTrans)
    if (!y %in% colnames(df)) return(NULL)
    df[[yLabel]] <- transformValues(df[[y]], yTrans)

    return(list(data=df, xLabel=xLabel, yLabel=yLabel))
}

#' Interface to modify the style of the plot points
#'
#' @param ns Namespace function
#' @param id Character: identifier
#' @param description Character: display text for user
#' @param help Character: extra text to help the user
#' @param colour Character: default colour
#' @param size Integer: default size
#' @param alpha Numeric: default transparency value
#'
#' @importFrom shiny tagList h4 helpText sliderInput
#'
#' @return HTML elements
#' @keywords internal
plotPointsStyle <- function(ns, id, description, help=NULL, size=2,
                            colour="black", alpha=1.0) {
    id2 <- function(att) ns(paste0(id, att))

    tagList(
        h4(description),
        if (!is.null(help)) helpText(help),
        sliderInput(id2("Size"), "Size", min=1, max=10, step=1, value=size,
                    width="100%"),
        colourInputMod(id2("Colour"), "Colour", value=colour),
        sliderInput(id2("Alpha"), "Opacity", min=0, max=1, step=0.01,
                    value=alpha, width="100%")
    )
}

#' Plot distribution using a density plot
#'
#' The tooltip shows the median, variance, maximum, minimum and number of non-NA
#' samples of each data series (if \code{data} contains names or column names,
#' those will be used as sample names and also appear in the tooltip).
#'
#' @param data Numeric, data frame or matrix: gene expression data or
#' alternative splicing event quantification values (sample names are based on
#' their \code{names} or \code{colnames})
#' @param groups List of sample names or vector containing the group name per
#' \code{data} value (read Details); if \code{NULL} or a character vector of
#' length 1, \code{data} values are considered from the same group
#' @param rug Boolean: show rug plot?
#' @param vLine Boolean: plot vertical lines (including descriptive statistics
#' for each group)?
#' @inheritDotParams stats::density.default -x -na.rm
#' @param title Character: plot title
#' @param psi Boolean: are \code{data} composed of PSI values? If \code{NULL},
#'   \code{psi = TRUE} if all \code{data} values are between 0 and 1
#' @param rugLabels Boolean: plot sample names in the rug?
#' @param rugLabelsRotation Numeric: rotation (in degrees) of rug labels; this
#'   may present issues at different zoom levels and depending on the proximity
#'   of \code{data} values
#' @param legend Boolean: show legend?
#' @param valueLabel Character: label for the value (by default, either
#' \code{Inclusion levels} or \code{Gene expression})
#'
#' @details Argument \code{groups} can be either:
#' \itemize{
#' \item{a list of sample names, e.g.
#' \code{list("Group 1"=c("Sample A", "Sample B"), "Group 2"=c("Sample C")))}}
#' \item{a character vector with the same length as \code{data}, e.g.
#' \code{c("Sample A", "Sample C", "Sample B")}.}
#' }
#'
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_plotOptions hc_tooltip
#' JS
#' @importFrom stats median var density
#'
#' @family functions to perform and plot differential analyses
#' @return \code{highchart} object with density plot
#' @export
#'
#' @examples
#' data   <- sample(20, rep=TRUE)/20
#' groups <- paste("Group", c(rep("A", 10), rep("B", 10)))
#' names(data) <- paste("Sample", 1:20)
#' plotDistribution(data, groups)
#'
#' # Using colours
#' attr(groups, "Colour") <- c("Group A"="pink", "Group B"="orange")
#' plotDistribution(data, groups)
plotDistribution <- function(data, groups=NULL, rug=TRUE, vLine=TRUE, ...,
                             title=NULL, psi=NULL, rugLabels=FALSE,
                             rugLabelsRotation=0, legend=TRUE,
                             valueLabel=NULL) {
    if (is.null(psi)) {
        psi <- min(data, na.rm=TRUE) >= 0 && max(data, na.rm=TRUE) <= 1
    }

    if (psi) {
        xMin   <- 0
        xMax   <- 1
        xLabel <- "Distribution of PSI values"
        id     <- "Inclusion level: "
    } else {
        xMin   <- NULL
        xMax   <- NULL
        xLabel <- "Distribution of gene expression"
        id     <- "Gene expression: "
    }
    if (!is.null(valueLabel)) id <- paste0(valueLabel, ": ")

    # Include X-axis zoom and hide markers
    hc <- highchart() %>%
        hc_chart(zoomType="x") %>%
        hc_xAxis(min=xMin, max=xMax, title=list(text=xLabel))

    if (legend) {
        autohideYaxis <- paste("
            function() {
                function isInvisible(el) { return !el.visible; }
                var groups = this.chart.legend.allItems,
                    areAllSeriesHidden = groups.map(isInvisible).every(Boolean);
                if (!areAllSeriesHidden) {
                    this.chart.yAxis[0].options.gridLineWidth = 1;
                    return this.value;
                } else {
                    this.chart.yAxis[0].options.gridLineWidth = 0;
                }
            }")
        hc <- hc %>% hc_yAxis(labels=list(formatter=JS(autohideYaxis)))
    }

    hc <- hc %>%
        hc_legend(enabled=legend) %>%
        hc_plotOptions(series = list(fillOpacity=0.3,
                                     marker=list(enabled=FALSE))) %>%
        hc_tooltip(
            headerFormat=NULL,
            pointFormat=paste(
                "{point.tooltipLabel}", br(),
                span(style="color:{point.color}", "\u25CF "),
                id, "{point.x:.2f}", br(),
                tags$b("{series.name}"), br(),
                "Number of samples: {series.options.samples}", br(),
                "Median: {series.options.median}", br(),
                "Variance: {series.options.var}", br(),
                "Range: {series.options.min} - {series.options.max}")) %>%
        export_highcharts()

    if (!is.null(title)) hc <- hc %>% hc_title(text=unname(title))

    if (is.null(groups)) {
        ns <- groups <- "All samples"
    } else if (is.list(groups)) {
        ns <- names(groups)
    } else {
        ns <- groups
    }

    count <- 0
    plotLines <- list()
    for (group in unique(ns)) {
        if (is.list(groups)) {
            filter <- groups[[group]]
        } else {
            filter <- groups == group
        }

        if (is.vector(data)) {
            row      <- data[filter]
        } else if (isTRUE(filter)) {
            row      <- data
        } else {
            filter   <- filter[filter %in% colnames(data)]
            row      <- data[ , filter]
        }

        # Prepare labels based on sample names (or the values themselves)
        if (!is.null(names(row))) {
            rugLabel     <- names(row)
            tooltipLabel <- rugLabel
        } else if (!is.null(colnames(row))) {
            rugLabel     <- colnames(row)
            tooltipLabel <- rugLabel
        } else {
            rugLabel     <- round(row, 2)
            tooltipLabel <- NULL
        }
        if (!rugLabels) rugLabel <- NULL

        row <- as.numeric(row)
        if (length(row) == 0) next

        # Stats
        med  <- roundDigits(median(row, na.rm=TRUE))
        vari <- roundDigits(var(row, na.rm=TRUE))
        max  <- roundDigits(max(row, na.rm=TRUE))
        min  <- roundDigits(min(row, na.rm=TRUE))
        samples <- sum(!is.na(row))

        colour <- unname(attr(groups, "Colour")[group])
        if (is.null(colour)) {
            colour <- JS(paste0("Highcharts.getOptions().colors[", count, "]"))
        }

        # Calculate the density of inclusion levels for each sample group
        den <- tryCatch(density(row, na.rm=TRUE, ...), error=return,
                        warning=return)
        if (length(row) == 1) {
            den  <- row
            vari <- max <- min <- med
        } else if (is(den, "error") || is(den, "warning")) {
            den <- NULL
        }

        if (!is.null(den)) {
            hc <- hc %>% hc_add_series(den, type="area", name=group, median=med,
                                       var=vari, samples=samples, max=max,
                                       color=colour, min=min)
            if (length(row) == 1) {
                len <- length(hc$x$hc_opts$series)
                hc$x$hc_opts$series[[len]]$visible <- FALSE
                if (legend) {
                    hc$x$hc_opts$series[[len]]$events$legendItemClick <- JS(
                        "function(e) { e.preventDefault() }")
                }
            }
        }
        # Rug plot
        if (rug) {
            isHexColour <- function(string) {
                # Explicitely ignores HEX colour codes with opacity
                grepl("^#{0,1}([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", string)
            }

            convertOpacityToHex <- function(opacity) {
                sprintf("%02x", round(opacity/100 * 255))
            }
            opacity <- convertOpacityToHex(60) # Opacity in percentage

            if (is(colour, "JS_EVAL")) {
                fill <- JS(sprintf("%s + \"%s\"", colour, opacity))
            } else if (isHexColour(colour)) {
                fill <- sprintf("%s%s", colour, opacity)
            } else {
                fill <- colour
            }

            # Add different, arbitrary y values per group (useful when only
            # displaying the rug plot)
            maxi <- ifelse("y" %in% names(den), max(den$y), max(den))
            y <- match(group, unique(ns)) * maxi / 1000
            hc <- hc %>%
                hc_scatter(
                    row, rep(y, length(row)), name=group, color=fill,
                    rugLabel=rugLabel, tooltipLabel=tooltipLabel,
                    marker=list(enabled=TRUE, radius=4, fillColor=fill),
                    median=med, var=vari, samples=samples, max=max, min=min)
            if (length(row) == 1) {
                len <- length(hc$x$hc_opts$series)
                hc$x$hc_opts$series[[len]] <- c(
                    hc$x$hc_opts$series[[len]],
                    median=med, var=vari, samples=samples, max=max, min=min)
            }
        }
        # Plot line with basic statistics
        if (vLine) {
            plotLines[[count + 1]] <- list(
                label=list(text=paste("Median:", med, "/ Variance:", vari)),
                color=colour, dashStyle="shortdash", width=2, value=med,
                zIndex=7)
        }
        count <- count + 1
    }

    # Add plotLines with information
    if (vLine) hc <- hc %>% hc_xAxis(plotLines = plotLines)

    # Show or hide rug labels
    rugSeries <- which(sapply(hc$x$hc_opts$series, "[[", "type") == "scatter")
    for (k in rugSeries) {
        hc$x$hc_opts$series[[k]]$dataLabels <- list(
            enabled=rugLabels,
            format="{point.rugLabel}",
            rotation=rugLabelsRotation)

        if (rugLabelsRotation != 0) {
            rugLabelsRotation <- rugLabelsRotation %% 360

            hc$x$hc_opts$series[[k]]$dataLabels$crop  <- FALSE
            if (rugLabelsRotation > 0 && rugLabelsRotation < 180) {
                align <- "right"
            } else if (rugLabelsRotation > 180 && rugLabelsRotation < 360) {
                align <- "left"
            } else {
                align <- "center"
            }
            hc$x$hc_opts$series[[k]]$dataLabels$align <- align
        }
    }
    return(hc)
}

#' Render boxplot
#'
#' @param data Data frame or matrix
#' @param outliers Boolean: draw outliers?
#' @param sortByMedian Boolean: sort box plots based on ascending median?
#' @param showXlabels Boolean: show labels in X axis?
#'
#' @importFrom reshape2 melt
#' @importFrom highcharter data_to_boxplot hc_add_series_list
#'
#' @return Box plot
#' @keywords internal
#'
#' @examples
#' psichomics:::renderBoxplot(data.frame(a=1:10, b=10:19, c=45:54))
renderBoxplot <- function(data, outliers=FALSE, sortByMedian=TRUE,
                          showXlabels=TRUE, title=NULL,
                          seriesName="Gene expression") {
    if (sortByMedian) {
        medians <- customColMedians(data, fast=TRUE)
        data <- data[ , order(medians)]
    }

    # Remove matrix rownames from melted data
    melted <- suppressMessages(melt(data))
    if (ncol(melted) == 3) {
        melted[[1]] <- NULL
        colnames(melted)[[1]] <- "variable"
    }
    value <- variable <- NULL
    dat <- data_to_boxplot(melted, value, variable, add_outliers=outliers)
    hc <- highchart() %>%
        hc_add_series_list(dat) %>%
        hc_chart(zoomType="x", type="column") %>%
        hc_plotOptions(boxplot=list(color="gray", fillColor="orange")) %>%
        hc_xAxis(type="category",
                 labels=list(enabled=showXlabels), visible=showXlabels) %>%
        hc_title(text=unname(title))
    if (min(melted$value) >= 0) hc <- hc %>% hc_yAxis(min=0)

    hc <- hc %>% export_highcharts()
    hc$x$hc_opts$series[[1]]$name <- seriesName
    return(hc)
}

#' Levene's test
#'
#' Performs a Levene's test to assess the equality of variances
#'
#' @details The implementation of this function is based on
#' \code{car:::leveneTest.default} with a more standard result.
#'
#' @param x Numeric vector or list of numeric vectors: non-numeric elements of a
#' list will be coerced with a warning
#' @param g Vector or factor: groups of elements in \code{x} (ignored with a
#' warning if \code{x} is a list)
#' @param centers Function used to calculate how much values spread; for
#' instance, \code{median} (default) or \code{mean}
#'
#' @importFrom stats complete.cases anova median lm
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic with a name describing it.}
#' \item{p.value}{the p-value for the test.}
#' \item{method}{the type of test applied.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @keywords internal
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

#' Create density sparklines for inclusion levels
#'
#' @inherit createSparklines
#' @param areSplicingEvents Boolean: are these splicing events (TRUE) or gene
#' expression (FALSE)?
#'
#' @importFrom highcharter highchart hc_credits hc_tooltip hc_chart hc_title
#' hc_xAxis hc_yAxis hc_exporting hc_legend hc_plotOptions
#' @importFrom shiny tags
#'
#' @keywords internal
createDensitySparklines <- function(data, events, areSplicingEvents=TRUE,
                                    groups=NULL, geneExpr=NULL,
                                    inputID="sparklineInput") {
    if (areSplicingEvents) {
        minX <- 0
        maxX <- 1
        id   <- "Inclusion levels"
        FUN  <- "showDiffSplicing"
    } else {
        minX <- NULL
        maxX <- NULL
        id   <- "Gene expression"
        FUN  <- "showDiffExpression"
    }

    hc <- highchart() %>%
        hc_tooltip(
            hideDelay=0, shared=TRUE, valueDecimals=getPrecision(),
            headerFormat=paste0(tags$small(paste0(id, ": {point.x:.2f}")),
                                tags$br()),
            pointFormat=paste(span(style="color:{point.color}", "\u25CF "),
                              tags$b("{series.name}"), br())) %>%
        hc_chart(width=120, height=20, backgroundColor="", type="areaspline",
                 margin=c(2, 0, 2, 0), style=list(overflow='visible')) %>%
        hc_title(text="") %>%
        hc_xAxis(visible=FALSE) %>%
        hc_credits(enabled=FALSE) %>%
        hc_yAxis(endOnTick=FALSE, startOnTick=FALSE, visible=FALSE) %>%
        hc_exporting(enabled=FALSE) %>%
        hc_legend(enabled=FALSE) %>%
        hc_plotOptions(series=list(cursor="non", animation=FALSE, lineWidth=1,
                                   marker=list(radius=1), fillOpacity=0.25))

    if (!is.null(minX) || !is.null(maxX))
        hc <- hc %>% hc_xAxis(min=minX, max=maxX)
    createSparklines(hc, data, events=events, groups=groups, geneExpr=geneExpr,
                     inputID=inputID)
}

#' Create sparkline charts to be used in a data table
#'
#' @param hc \code{highchart} object
#' @param data Character: HTML-formatted data series of interest
#' @param id Character: Shiny input identifier
#' @param events Character: event identifiers
#' @param groups Character: name of the groups used for differential analyses
#' @param geneExpr Character: name of the gene expression dataset
#' @param inputID Character: identifier of input to get attributes of clicked
#'   event (Shiny only)
#'
#' @importFrom jsonlite toJSON
#'
#' @return HTML element with sparkline data
#' @keywords internal
createSparklines <- function(hc, data, events, groups=NULL, geneExpr=NULL,
                             inputID="sparklineInput", ...) {
    hc <- as.character(toJSON(hc$x$hc_opts, auto_unbox=TRUE))
    hc <- substr(hc, 1, nchar(hc)-1)

    if (!is.null(groups) && !identical(groups, "")) {
        groups <- toJSarray(unlist(groups))
    } else {
        groups <- "null"
    }
    if (is.null(geneExpr) || geneExpr == "") {
        geneExpr <- ""
        type     <- "event"
    } else {
        geneExpr <- sprintf(", geneExpr: '%s'", geneExpr)
        type     <- "gene"
    }
    params  <- sprintf("{%s: '%s', groups: %s%s}",
                       type, events, groups, geneExpr)
    onclick <- sprintf("Shiny.setInputValue('%s', %s, {priority: 'event'});",
                       inputID, params)
    json <- paste0(hc, ',"series":', data, "}")
    sparklines <- sprintf(paste(
        '<sparkline onclick="%s"',
        'style="cursor:pointer;" data-sparkline=\'%s\'/>'),
        onclick, json)
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
#' @param threshold Integer: minimum number of values per group
#' @param analyses Character: analyses to perform (see Details)
#' @param step Numeric: number of events before the progress bar is updated
#' (a bigger number allows for a faster execution)
#'
#' @importFrom stats kruskal.test median wilcox.test t.test var density
#' @importFrom methods is
#'
#' @return A row from a data frame with the results
#' @keywords internal
singleDiffAnalyses <- function(vector, group, threshold=1, step=100,
                               analyses=c("wilcoxRankSum", "ttest", "kruskal",
                                          "levene", "fligner")) {
    colour  <- attr(group, "Colour")
    series  <- split(vector, group)
    samples <- vapply(series, function(i) sum(!is.na(i)), integer(1))
    valid   <- names(series)[samples >= threshold]
    if(length(valid) == 0) return(NULL)

    inGroup <- group %in% valid
    group   <- group[inGroup]
    vector  <- vector[inGroup]
    len     <- length(valid)

    # Variance and median
    med <- lapply(series, median, na.rm=TRUE)
    var <- lapply(series, var, na.rm=TRUE)

    # Improve display of group names
    names(samples) <- paste0("(", names(samples), ")")
    names(med) <- paste0("(", names(med), ")")
    names(var) <- paste0("(", names(var), ")")

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
        wilcox  <- suppressWarnings(
            tryCatch(wilcox.test(vector[typeOne], vector[!typeOne]),
                     error=return))
        if (is(wilcox, "error")) wilcox <- NULL
        # } else if (any("wilcoxSignedRank" == analyses) && len == 1) {
        #     # Wilcoxon signed rank test
        #     wilcox <- suppressWarnings(wilcox.test(vector))
    }

    # Kruskal-Wallis test (2 or more groups)
    kruskal <- NULL
    if (any("kruskal" == analyses) && len >= 2) {
        kruskal <- tryCatch(kruskal.test(vector, group), error=return)
        if (is(kruskal, "error")) kruskal <- NULL
    }

    # Levene's test (2 or more groups)
    levene <- NULL
    if (any("levene" == analyses) && len >= 2) {
        levene <- suppressWarnings(
            tryCatch(leveneTest(vector, group), error=return))
        if (is(levene, "error")) levene <- NULL
    }

    # Fligner-Killeen test (2 or more groups)
    fligner <- NULL
    if (any("fligner" == analyses) && len >= 2) {
        fligner <- suppressWarnings(
            tryCatch(fligner.test(vector, group), error=return))
        if (is(fligner, "error") || is.infinite(fligner$statistic))
            fligner <- NULL
    }

    # Density sparklines
    sparkline <- NULL
    if (any("density" == analyses)) {
        data         <- NULL
        validSeries  <- series[valid]
        groupsColour <- NULL
        for (each in seq(validSeries)) {
            group <- validSeries[[each]]
            # Calculate data density for each sample group with a greatly
            # reduced number of points for faster execution
            den <- tryCatch(density(group, n=10, na.rm=TRUE), error=return,
                            warning=return)
            if (is(den, "error") || is(den, "warning"))
                data <- c(data, "")
            else
                data <- c(data, paste(sprintf('{"x":%s,"y":%s}', den$x, den$y),
                                      collapse=","))

            if (!is.null(colour))
                groupsColour <- c(groupsColour,
                                  unname(colour[names(validSeries)[[each]]]))
        }
        if (!is.null(colour)) {
            sparkline <- paste(
                sprintf('{"name":"%s", "data":[%s], "color":"%s"}',
                        names(validSeries), data, groupsColour), collapse=",")
        } else {
            sparkline <- paste(
                sprintf('{"name":"%s", "data":[%s]}',
                        names(validSeries), data), collapse=",")
        }
        sparkline <- paste("[", sparkline, "]")
    }

    vector <- c("Distribution"=sparkline,
                "Survival by PSI cutoff"=as.numeric(NA),
                "Optimal PSI cutoff"=as.numeric(NA),
                "Log-rank p-value"=as.numeric(NA),
                "Samples"=samples, "T-test"=ttest,
                "Wilcoxon"=wilcox, "Kruskal"=kruskal, "Levene"=levene,
                "Fligner-Killeen"=fligner, "Variance"=var, "Median"=med)
    vector <- vector[!vapply(vector, is.null, logical(1))] # Remove NULL
    return(vector)
}

#' @importFrom fastmatch fmatch
#' @importFrom plyr rbind.fill
convertListOfLists2DataFrame <- function(stats) {
    # Check the column names of the different columns
    ns    <- lapply(stats, names)
    uniq  <- unique(ns)
    match <- fmatch(ns, uniq)

    ll <- lapply(stats, function(i) lapply(i, unname))
    ll <- lapply(ll, unlist)
    ldf <- lapply(seq_along(uniq), function(k) {
        elems <- match == k
        df2   <- t(data.frame(ll[elems], stringsAsFactors=FALSE))
        cols  <- colnames(df2)
        if (nrow(df2) == 0) return(NULL)

        df2           <- data.frame(df2, stringsAsFactors=FALSE)
        colnames(df2) <- cols
        rownames(df2) <- names(stats)[elems]
        return(df2)
    })
    df <- rbind.fill(ldf)
    rownames(df) <- unlist(lapply(ldf, rownames))
    return(df)
}

convertNumberCols2Numbers <- function(df) {
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
        int <- apply(df[ , numericCols, drop=FALSE], 2,
                     function(i) all(is.whole(i), na.rm=TRUE))
        intCols <- numericCols
        intCols[numericCols] <- int
        if (any(intCols)) {
            df[ , intCols] <- apply(df[ , intCols, drop=FALSE], 2, as.integer)
        }
    }
    return(df)
}

combineSplicingEventInfo <- function(df, data) {
    events   <- rownames(df)
    info     <- parseSplicingEvent(events, data=data, pretty=TRUE)
    if (is.null(info)) return(df)
    infoGene <- prepareGenePresentation(info$gene)
    df       <- cbind("Event type"=info$subtype, "Chromosome"=info$chrom,
                      "Strand"=info$strand, "Gene"=unlist(infoGene), df)
    rownames(df) <- events
    return(df)
}

#' Perform statistical analyses
#'
#' @param data Data frame or matrix: gene expression or alternative splicing
#' quantification
#' @param groups Named list of characters (containing elements belonging to each
#' group) or character vector (containing the group of each individual sample);
#' if \code{NULL}, sample types are used instead when available, e.g. normal,
#' tumour and metastasis
#' @param analyses Character: statistical tests to perform (see Details)
#' @param pvalueAdjust Character: method used to adjust p-values (see Details)
#' @param geneExpr Character: name of the gene expression dataset (only required
#' for density sparklines available in the interactive mode)
#' @inherit createDensitySparklines
#'
#' @importFrom stats p.adjust
#'
#' @details
#' The following statistical analyses may be performed simultaneously via the
#' \code{analysis} argument:
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
#' The following p-value adjustment methods are supported via the
#' \code{pvalueAdjust} argument:
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
#' @family functions to perform and plot differential analyses
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
diffAnalyses <- function(data, groups=NULL,
                         analyses=c("wilcoxRankSum", "ttest", "kruskal",
                                    "levene", "fligner"),
                         pvalueAdjust="BH", geneExpr=NULL,
                         inputID="sparklineInput") {
    # cl <- parallel::makeCluster(getOption("cl.cores", getCores()))
    step <- 50 # Avoid updating progress too frequently
    updateProgress("Performing statistical analysis",
                   divisions=5 + round(nrow(data)/step))
    time <- Sys.time()

    if (is.null(groups)) {
        ids    <- names(data)
        groups <- parseTCGAsampleTypes(ids)
    } else if (is.list(groups)) {
        groups <- discardOutsideSamplesFromGroups(groups, colnames(data))
        data   <- data[ , unlist(groups)]

        colour <- attr(groups, "Colour")
        groups <- rep(names(groups), sapply(groups, length))
        attr(groups, "Colour") <- colour
    }
    originalGroups <- unique(groups)
    if (identical(originalGroups, "All samples")) originalGroups <- NULL

    # Prepare groups with respective colours
    colour <- attr(groups, "Colour")
    groups <- factor(groups)
    if ( !is.null(colour) ) attr(groups, "Colour") <- colour

    count <- 0
    stats <- apply(data, 1, function(...) {
        count <<- count + 1
        if (count %% step == 0)
            updateProgress("Performing statistical analysis", console=FALSE)
        return(singleDiffAnalyses(...))
    }, groups, threshold=1, step=step, analyses=analyses)
    display(Sys.time() - time)

    updateProgress("Preparing data")
    time <- Sys.time()
    df   <- convertListOfLists2DataFrame(stats)
    df   <- convertNumberCols2Numbers(df)
    display(Sys.time() - time)

    # Calculate delta variance and delta median if there are only 2 groups
    deltaVar <- df[ , grepl("Variance", colnames(df)), drop=FALSE]
    if (ncol(deltaVar) == 2) {
        updateProgress("Calculating delta variance and median")
        time     <- Sys.time()
        deltaVar <- deltaVar[ , 1] - deltaVar[ , 2]
        deltaMed <- df[, grepl("Median", colnames(df))]
        deltaMed <- deltaMed[ , 1] - deltaMed[ , 2]
        # Avoid R warnings in Windows by setting Unicode codes in column names
        colns <- colnames(df)
        df    <- cbind(df, deltaVar, deltaMed)
        colnames(df) <- c(colns, "\u2206 Variance", "\u2206 Median")
        display(Sys.time() - time)
    }

    if (any(pvalueAdjust == c("BH", "BY", "bonferroni", "holm", "hochberg",
                              "hommel"))) {
        updateProgress("Adjusting p-values", detail=pvalueAdjust)
        cols <- grep("p.value", colnames(df))[-1]
        if (length(cols > 0)) {
            time <- Sys.time()
            pvalue <- df[cols]
            adjust <- apply(pvalue, 2, p.adjust, pvalueAdjust)
            names  <- paste0(colnames(pvalue), " (", pvalueAdjust, " adjusted)")

            if (!is.matrix(adjust)) adjust <- t(as.matrix(adjust))
            colnames(adjust) <- names

            # Place the adjusted p-values next to the respective p-values
            len <- ncol(df)
            order <- seq(len)
            for (i in seq_along(cols))
                order <- append(order, len+i, after=which(order == cols[i]))
            df <- cbind(df, adjust)[order]
            display(Sys.time() - time)
        }
    }

    areEvents <- areSplicingEvents(rownames(df), data=data)
    if (areEvents) {
        updateProgress("Including splicing event information")
        df <- combineSplicingEventInfo(df, data)
    }

    if (any("density" == analyses)) {
        updateProgress("Preparing density plots")
        time <- Sys.time()

        df[ , "Distribution"] <- createDensitySparklines(
            df[ , "Distribution"], rownames(df), areEvents,
            groups=originalGroups, geneExpr=geneExpr, inputID=inputID)
        name <- ifelse(areEvents, "PSI.distribution", "GE.distribution")
        colnames(df)[match("Distribution", colnames(df))] <- name
        display(Sys.time() - time)
    }

    # Properly set column names
    col <- colnames(df)
    col <- gsub(".", " ", col, fixed=TRUE)
    col <- gsub("p value", "p-value", col, fixed=TRUE)
    colnames(df) <- col

    # parallel::stopCluster(cl)
    attr(df, "rowData") <- getSplicingEventData(data)
    df <- preserveAttributes(df)
    return(df)
}

prettifyEventID <- function(event, data=NULL) {
    eventData <- findEventData(event, data=data)
    hasID     <- !is.null(eventData$id)
    parsed    <- parseSplicingEvent(event, char=!hasID, data=data)
    if (is.null(parsed)) {
        parsed <- event
    } else if (hasID) {
        parsed <- parsed$id
    }
    return(parsed)
}

#' Set of functions to render differential analyses (plot and table)
#'
#' @inheritParams appServer
#' @param analysesType Character: type of analyses (\code{GE} or \code{PSI})
#' @param analysesID Character: identifier
#' @param getAnalysesData Function: get analyses data
#' @param getAnalysesFiltered Function: get filtered analyses data
#' @param setAnalysesFiltered Function: set filtered analyses data
#' @param getAnalysesSurvival Function: get survival data
#' @param getAnalysesColumns Function: get columns
#' @param setAnalysesColumns Function: set columns
#' @param getResetPaging Function: get toggle of reset paging
#' @param setResetPaging Function: set toggle of reset paging
#'
#' @importFrom DT dataTableProxy selectRows replaceData
#' @importFrom shinyjs toggleElement toggleState
#' @importFrom utils write.table
#'
#' @inherit psichomics return
#' @keywords internal
analysesTableSet <- function(session, input, output, analysesType, analysesID,
                             getAnalysesData, getAnalysesFiltered,
                             setAnalysesFiltered, getAnalysesSurvival,
                             getAnalysesColumns, setAnalysesColumns,
                             getResetPaging, setResetPaging) {
    # Save selected points in the table
    observe({
        selected <- input$statsTable_rows_selected
        setSelectedPoints(analysesID, selected)
    })

    if (analysesType == "PSI") {
        searchableCols <- 5
        visibleCols    <- 6:8
    } else if (analysesType == "GE") {
        searchableCols <- 1
        visibleCols    <- 5:6
    }

    # Render table with sparklines
    output$statsTable <- renderDataTableSparklines({
        stats <- getAnalysesFiltered()
        if (!is.null(stats)) {
            # Discard columns of no interest
            cols <- colnames(stats)
            cols <- cols[!grepl("method|data.name", cols)]
            setAnalysesColumns(cols)
            return(stats[ , cols])
        }
    }, style="bootstrap", filter="top", server=TRUE, extensions="Buttons",
    options=list(pageLength=10, dom="Bfrtip", buttons=I("colvis"),
                 columnDefs=list(
                     list(targets=searchableCols, searchable=FALSE),
                     list(targets=visibleCols, visible=FALSE))))

    # Update table with filtered information
    proxy <- dataTableProxy("statsTable")
    observe({
        stats <- getAnalysesData()
        if (!is.null(stats)) {
            # Bind preview of survival curves based on value cutoff
            optimSurv <- getAnalysesSurvival()
            if (!is.null(optimSurv)) {
                survCols <- sprintf(c("Optimal %s cutoff", "Log-rank p-value",
                                      "Survival by %s cutoff"), analysesType)
                stats[[survCols[[1]]]] <- optimSurv[[1]]
                stats[[survCols[[2]]]] <- optimSurv[[2]]
                stats[[survCols[[3]]]] <- optimSurv[[3]]
            }

            # Filter by highlighted events and events in the zoomed area
            events  <- getHighlightedPoints(analysesID)
            zoom    <- getZoom(analysesID)

            zoomed <- NULL
            if (!is.null(zoom)) {
                x <- input$xAxis
                y <- input$yAxis
                if (!is.null(x) && !is.null(y)) {
                    res <- transformData(input, stats, x, y)
                    if (!is.null(res)) {
                        stats  <- res$data
                        xLabel <- res$xLabel
                        yLabel <- res$yLabel

                        xStats <- stats[[xLabel]]
                        xZoom  <- zoom$xmin <= xStats & xStats <= zoom$xmax
                        yStats <- stats[[yLabel]]
                        yZoom  <- zoom$ymin <= yStats & yStats <= zoom$ymax
                        zoomed <- intersect(which(xZoom), which(yZoom))
                    }
                }
            }

            # Filter rows based on highlighted and/or zoomed in events
            if (!is.null(events) && !is.null(zoomed)) {
                rowFilter <- intersect(events, zoomed)
            } else if (!is.null(events)) {
                rowFilter <- events
            } else if (!is.null(zoomed)) {
                rowFilter <- zoomed
            } else {
                rowFilter <- TRUE
            }
            stats <- stats[rowFilter, ]

            # Keep previously selected rows if possible
            before   <- isolate(getAnalysesFiltered())
            selected <- isolate(input$statsTable_rows_selected)
            selected <- rownames(before)[isolate(selected)]
            selected <- which(rownames(stats) %in% selected)
            if (length(selected) < 1) selected <- NULL

            # Set new data
            setAnalysesFiltered(stats)

            # Keep cols or no data will be rendered
            cols  <- getAnalysesColumns()
            cols  <- cols[cols %in% colnames(stats)]
            stats <- stats[ , cols]

            # Check if paging should be reset
            resetPaging <- isolate(getResetPaging())
            if (is.null(resetPaging)) resetPaging <- TRUE
            setResetPaging(resetPaging)

            # Round numbers based on significant digits
            cols <- colnames(stats)
            type <- sapply(cols, function(i) class(stats[[i]]))
            numericCols <- cols[type == "numeric"]

            # Round numbers based on significant digits
            if (nrow(stats) > 0) {
                for (col in numericCols) {
                    stats[ , col] <- suppressWarnings(
                        as.numeric(signifDigits(stats[ , col])))
                }
            }
            replaceData(proxy, stats, resetPaging=resetPaging,
                        clearSelection="none")
        }
    })

    # Hide table toolbar if statistical table is not displayed
    observe(toggleElement(
        "tableToolbar", condition=!is.null(getAnalysesData())))

    # Discard columns from data frame containing information to render plots
    discardPlotsFromTable <- function(df) {
        plotCols <- TRUE
        if (!is.null(df)) {
            ns <- sprintf(c("%s distribution", "Survival by %s cutoff"),
                          analysesType)
            plotCols <- -match(ns, colnames(df))
            plotCols <- plotCols[!is.na(plotCols)]
            if (length(plotCols) == 0) plotCols <- TRUE
        }
        return(df[ , plotCols])
    }

    if (analysesType == "PSI") {
        filenameText <- "Differential splicing analyses"
        rownamesCol  <- "AS event"
    } else if (analysesType == "GE") {
        filenameText <- "Differential expression analyses"
        rownamesCol  <- "Gene"
    }

    # Download whole table
    output$downloadAll <- downloadHandler(
        filename=paste(getCategory(), filenameText),
        content=function(file) {
            stats <- getAnalysesData()
            stats <- discardPlotsFromTable(stats)

            # Include updated survival analyses
            optimSurv <- getAnalysesSurvival()
            if (!is.null(optimSurv)) {
                cols <- sprintf(c("Optimal %s cutoff", "Log-rank p-value"),
                                analysesType)
                stats[[cols[[1]]]] <- optimSurv[[1]]
                stats[[cols[[2]]]] <- optimSurv[[2]]
            }

            stats <- cbind(rownames(stats), stats)
            colnames(stats)[[1]] <- rownamesCol
            write.table(stats, file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )

    # Download filtered table
    output$downloadSubset <- downloadHandler(
        filename=paste(getCategory(), filenameText),
        content=function(file) {
            stats <- getAnalysesFiltered()
            stats <- discardPlotsFromTable(stats)
            stats <- stats[input$statsTable_rows_all, ]

            stats <- cbind(rownames(stats), stats)
            colnames(stats)[[1]] <- rownamesCol
            write.table(stats, file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )

    # Create groups based on a given filter
    groupBasedOnAnalysis <- function(filter, description="") {
        stats <- getAnalysesFiltered()
        stats <- discardPlotsFromTable(stats)
        stats <- stats[filter, ]

        if (analysesType == "PSI") {
            ASevents <- rownames(stats)
            genes  <- unique(getGenesFromSplicingEvents(ASevents, data=stats))
            origin <- "Selection from differential splicing analysis"
            title  <- "Differential splicing selection"
        } else if (analysesType == "GE") {
            genes <- rownames(stats)

            ASevents <- getASevents()
            if ( !is.null(ASevents) ) {
                ASevents <- getSplicingEventFromGenes(genes, ASevents)
            } else {
                ASevents <- character(0)
            }
            origin <- "Selection from differential expression analysis"
            title  <- "Differential expression selection"
        }

        group <- cbind("Names"=title, "Subset"=origin, "Input"=origin,
                       "ASevents"=list(ASevents), "Genes"=list(genes))
        appendNewGroups("ASevents", group)
        infoModal(
            session, title="New group created",
            "New group created", description, "and containing:",
            div(style="font-size: 22px;", length(ASevents), "splicing events"),
            div(style="font-size: 22px;", length(genes), "genes"))
    }

    if (analysesType == "PSI")     groupedElements <- "splicing events"
    else if (analysesType == "GE") groupedElements <- "genes"
    groupsText <- sprintf(c("based on the %s shown in the table",
                            "based on selected %s"), groupedElements)

    # Create groups based on splicing events displayed in the table
    observeEvent(input$groupByDisplayed, groupBasedOnAnalysis(
        input$statsTable_rows_all, groupsText[[1]]))

    # Create groups based on selected splicing events
    observeEvent(input$groupBySelected, groupBasedOnAnalysis(
        input$statsTable_rows_selected, groupsText[[2]]))

    # Disable groups based on selected AS events when no groups are selected
    observe(toggleState("groupBySelectedContainer",
                        !is.null(input$statsTable_rows_selected)))
}

#' Set up environment and redirect user to a page based on click information
#'
#' @param click List: click information
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param survival Boolean: redirect to survival page?
#'
#' @rdname analysesTableSet
#'
#' @keywords internal
processClickRedirection <- function(click, psi=NULL, survival=FALSE) {
    if (is.null(click)) return(NULL)

    groups <- ifelse(is.null(click$groups),
                     "null", toJSarray(click$groups))
    if (groups == "['']") groups <- "null"

    isSplicingEvent <- !is.null(click$event)
    if (isSplicingEvent) setASevent(click$event, data=psi)

    if (!survival) {
        if (isSplicingEvent) {
            js <- sprintf("showDiffSplicing(%s);", groups)
        } else {
            js <- sprintf("showDiffExpression('%s', %s, '%s');",
                          click$gene, groups, click$geneExpr)
        }
    } else {
        js <- sprintf("showSurvCutoff('%s', %s, autoParams = true, psi = %s)",
                      click[[1]], groups,
                      ifelse(isSplicingEvent, "true", "false"))
    }
    runjs(js)
}

#' @rdname analysesTableSet
#'
#' @importFrom stringr str_split
#' @importFrom shinyjs toggleState
analysesPlotSet <- function(session, input, output, analysesType, analysesID,
                            getAnalysesData, getAnalysesFiltered,
                            getAnalysesSurvival) {
    ns <- session$ns

    # Toggle visibility of elements regarding event options
    observe({
        stats <- getAnalysesData()
        if (is.null(stats)) {
            show("missingDiffAnalyses")
            hide("eventOptions")
        } else {
            hide("missingDiffAnalyses")
            show("eventOptions")
        }
    })

    # Update columns available to plot
    observe({
        stats <- getAnalysesData()
        optimSurv <- getAnalysesSurvival()
        if (!is.null(optimSurv)) {
            optimSurvCols <- sprintf(
                c("Optimal %s cutoff", "Log rank p-value"), analysesType)
            stats[[optimSurvCols[1]]] <- optimSurv[[1]]
            stats[[optimSurvCols[2]]] <- optimSurv[[2]]

            # Show these columns at the end
            names <- colnames(stats)
            colsMatch <- match(optimSurvCols, names)
            stats <- stats[c(names[-colsMatch], names[colsMatch])]
        }
        eventPlotOptions(session, stats, isolate(input$xAxis),
                         isolate(input$yAxis), isolate(input$labelSortBy))
    })

    # Update alternative splicing events and genes available to label
    observe({
        diffAnalyses <- getAnalysesData()
        if (!is.null(diffAnalyses)) {
            if (analysesType == "PSI") {
                ASevents <- rownames(diffAnalyses)
                # TODO: use prettifyEventID()?
                names(ASevents) <- gsub("_", " ", ASevents)
                updateSelectizeInput(session, "labelEvents", server=TRUE,
                                     choices=ASevents, selected=character(0))
                allGenes <- sort(unique(unlist(
                    str_split(diffAnalyses$Gene, "/"))))
                updateSelectizeInput(session, "labelGenes", server=TRUE,
                                     choices=allGenes, selected=character(0))
            } else if (analysesType == "GE") {
                updateSelectizeInput(session, "labelGenes", server=TRUE,
                                     choices=rownames(diffAnalyses),
                                     selected=character(0))
            }
        } else {
            if (analysesType == "PSI") {
                updateSelectizeInput(session, "labelEvents", server=TRUE,
                                     choices=character(0),
                                     selected=character(0))
            }

            updateSelectizeInput(session, "labelGenes", server=TRUE,
                                 choices=character(0), selected=character(0))
        }
    })

    # Interface elements to highlight values in the plot
    lapply(c("x", "y"), function(axis) {
        observe({
            highlightUI <- function(label, min, max) {
                highlightId <- ns(paste0(label, "Highlight"))
                sliderMinId <- ns(paste0(label, "SliderMin"))
                sliderMaxId <- ns(paste0(label, "SliderMax"))
                sliderInvId <- ns(paste0(label, "SliderInv"))

                # Round max and min numbers with two decimal points
                max <- ceiling(max*100)/100
                min <- floor(min*100)/100

                conditionalPanel(
                    sprintf("input[id='%s']", highlightId),
                    fluidRow(
                        column(6, textInput(sliderMinId, "Lower limit \u2265",
                                            placeholder=min, width="100%")),
                        column(6, textInput(sliderMaxId, "Upper limit \u2264",
                                            placeholder=max, width="100%"))),
                    checkboxInput(sliderInvId, "Invert highlighted values"),
                    helpText("The data in the table is also filtered",
                             "according to highlighted events."))
            }

            stats <- getAnalysesData()
            optimSurv <- getAnalysesSurvival()
            if (!is.null(optimSurv)) {
                cols <- sprintf(c("Optimal %s cutoff", "Log-rank p-value"),
                                analysesType)
                stats[[cols[[1]]]] <- optimSurv[[1]]
                stats[[cols[[2]]]] <- optimSurv[[2]]
            }

            value <- input[[paste0(axis, "Axis")]]
            if (is.null(stats) || is.null(value)) {
                output[[paste0(axis, "HighlightValues")]] <- renderUI(NULL)
                return(NULL)
            }

            trans <- input[[paste0(axis, "Transform")]]
            label <- transformOptions(value, trans)
            if (!value %in% colnames(stats)) {
                output[[paste0(axis, "HighlightValues")]] <- renderUI(NULL)
                return(NULL)
            }
            vals  <- transformValues(stats[[value]], trans)
            rangeNos <- range(vals, na.rm=TRUE)
            minNo    <- min(rangeNos)
            maxNo    <- max(rangeNos)

            output[[paste0(axis, "HighlightValues")]] <- renderUI(
                highlightUI(axis, minNo, maxNo) )
        })
    })

    # Disable labelling elements as appropriate
    observe(toggleState("labelTopOptions", input$labelTopEnable))
    observe(toggleState("labelGenes", input$labelGeneEnable))
    if (analysesType == "PSI")
        observe(toggleState("labelEvents", input$labelEventEnable))

    # Prepare labelled points
    observeEvent(input$labelPoints, {
        isolate({
            labelTopEnable   <- input$labelTopEnable
            labelEventEnable <- input$labelEventEnable # Only for PSIs
            labelGeneEnable  <- input$labelGeneEnable
            labelSortBy      <- input$labelSortBy
            labelTop         <- input$labelTop
            labelOrder       <- input$labelOrder == "decreasing"
            events           <- input$labelEvents # Only for PSIs
            genes            <- input$labelGenes
            displayGene      <- input$labelOnlyGene # Only for PSIs
            diffAnalyses     <- getAnalysesData()
        })

        labelled <- NULL
        # Label top events
        if (labelTopEnable && !identical(labelSortBy, "")) {
            sorted   <- order(diffAnalyses[[labelSortBy]],
                              decreasing=labelOrder)
            labelled <- head(sorted, labelTop)
        }

        if (analysesType == "PSI") {
            # Label selected alternative splicing events
            if (labelEventEnable && !identical(events, ""))
                labelled <- c(labelled, match(events, rownames(diffAnalyses)))

            # Label alternative splicing events based on selected genes
            if (labelGeneEnable && !identical(genes, "")) {
                # Unlist all genes and save original indexes to quickly find
                # genes across multiple alternative splicing events
                allGenes <- str_split(diffAnalyses$Gene, "/")
                len      <- sapply(allGenes, length)
                genes    <- sprintf("^%s$", genes)
                match    <- lapply(genes, grep, unlist(allGenes))
                index    <- rep(seq(diffAnalyses$Gene), len)[unlist(match)]
                labelled <- c(labelled, unique(index))
            }
        } else if (analysesType == "GE") {
            # Label selected genes
            if (labelGeneEnable && !identical(genes, ""))
                labelled <- c(labelled, match(genes, rownames(diffAnalyses)))
        }

        attr(labelled, "displayOnlyGene") <- displayGene
        setLabelledPoints(analysesID, labelled)
    })

    # Unlabel points
    observeEvent(input$unlabelPoints, setLabelledPoints(analysesID, NULL))

    # Plot data and prepare tooltip
    observe({
        stats     <- getAnalysesData()
        filtered  <- getAnalysesFiltered()
        eventData <- attr(stats, "eventData")

        x <- input$xAxis
        y <- input$yAxis
        if (is.null(stats) || is.null(x) || is.null(y)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }

        # Include survival data
        optimSurv <- getAnalysesSurvival()
        if (!is.null(optimSurv)) {
            cols <- sprintf(c("Optimal %s cutoff", "Log-rank p-value"),
                            analysesType)
            stats[[cols[[1]]]] <- optimSurv[[1]]
            stats[[cols[[2]]]] <- optimSurv[[2]]
        }

        res <- transformData(input, stats, x, y)
        if (is.null(res)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }
        statsDf <- res$data
        xLabel  <- res$xLabel
        yLabel  <- res$yLabel

        ggplotServer(
            input=input, output=output, id=analysesID, df=statsDf, x=xLabel,
            y=yLabel, eventData=stats, plot={
                parseHighlight <- function(input, arg) {
                    argStr <- function(...) paste0(arg, ...)

                    if (!input[[argStr("Highlight")]]) return(NULL)

                    highlightMin <- input[[argStr("SliderMin")]]
                    highlightMax <- input[[argStr("SliderMax")]]

                    # Parse string and expect a numeric value
                    evalString <- function(arg) {
                        value <- tryCatch(
                            suppressWarnings(eval(parse(text=arg),
                                                  envir=baseenv())),
                            error=return)
                        if (!is.numeric(value) || is(value, "error"))
                            value <- NULL
                        return(value)
                    }

                    highlightMax <- evalString(highlightMax)
                    highlightMin <- evalString(highlightMin)

                    noMin <- is.null(highlightMin)
                    noMax <- is.null(highlightMax)
                    minLTmax <- !noMin && !noMax &&
                        isTRUE(highlightMin >= highlightMax)
                    if ((noMin && noMax) || minLTmax) return(NULL)

                    highlight <- c(lower=highlightMin, upper=highlightMax)
                    attr(highlight, "inverted") <- input[[argStr("SliderInv")]]
                    return(highlight)
                }

                highlightX <- parseHighlight(input, "x")
                highlightY <- parseHighlight(input, "y")

                # Check selected events
                selected <- getSelectedPoints(analysesID)
                selected <- rownames(filtered)[selected]
                selected <- which(rownames(statsDf) %in% selected)
                if (length(selected) < 1) selected <- NULL

                params <- list(size=input$baseSize, col=input$baseColour,
                               alpha=input$baseAlpha)
                highlightParams <- list(size=input$highlightedSize,
                                        col=input$highlightedColour,
                                        alpha=input$highlightedAlpha)
                selectedParams  <- list(size=input$selectedSize,
                                        col=input$selectedColour,
                                        alpha=input$selectedAlpha)
                labelledParams  <- list(size=input$labelledSize,
                                        col=input$labelledColour,
                                        alpha=input$labelledAlpha)

                zoom <- getZoom(analysesID)
                if (!is.null(zoom)) {
                    xlim <- c(zoom$xmin, zoom$xmax)
                    ylim <- c(zoom$ymin, zoom$ymax)
                } else {
                    xlim <- NULL
                    ylim <- NULL
                }

                labelled <- getLabelledPoints(analysesID)
                eventPlot <- createEventPlotting(
                    statsDf, xLabel, yLabel, params, highlightX, highlightY,
                    highlightParams, selected, selectedParams, labelled,
                    labelledParams, xlim=xlim, ylim=ylim)
                setHighlightedPoints(analysesID, eventPlot$highlighted)
                eventPlot$plot[[1]]
            })
    })

    ggplotAuxServer(input, output, analysesID)
}

attr(analysesUI, "loader") <- "app"
attr(analysesServer, "loader") <- "app"
